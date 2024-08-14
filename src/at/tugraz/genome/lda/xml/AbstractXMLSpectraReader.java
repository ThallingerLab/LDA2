/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2021 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp 
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */

package at.tugraz.genome.lda.xml;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Vector;
import java.util.zip.GZIPInputStream;
import java.util.zip.Inflater;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;

import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.maspectras.quantification.CgDefines;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgScan;

/**
 * StAX based mzXML Reader for the Package.
 * 
 * This class reads an mzXML file with hierarchical scans. It uses two methods
 * defined in the AddScan interface to return created header and spectrum
 * information to the caller.
 * 
 * @author Leonida M. Lamp
 * @author Juergen Hartler
 */
public abstract class AbstractXMLSpectraReader implements XMLSpectraReader
{
  
  //the accepted file types
  public final static String FILE_TYPE_MZ_XML = "mzXML";
  public final static String FILE_TYPE_MZ_ML = "mzML";
	
  //the reader
  XMLStreamReader reader_;
  
  //the reader can serve several AddScan interfaces
  AddScan[] adders_;
  
  // if true MS/MS spectra are parsed too
  private boolean parseMsMs_;
  
  // the provided input stream
  private InputStream inStream_;
  
  // the multiplication factor to be used, to create integer values out of the m/z float values
  private int multiplicationFactorForInt_ = CgDefines.mzMultiplicationFactorForInt;
  // constant to give extra digits to the m/z value in integer format
  public static final int ONE_MILLION = 1000000;
  // the lowest available m/z value in integer format (original m/z value times the multiplication factor for int)
  private int lowestMz_ = ONE_MILLION * this.multiplicationFactorForInt_;
  // the highest available m/z value in integer format (original m/z value times the multiplication factor for int)
  private int highestMz_ = 0;
  
  //the maximally allowed m/z range provided by all the AddScan interfaces
  private Range maxRange_ = null;
  
  /** the polarity of the current scan*/
  protected int currentPolarity_ = CgDefines.POLARITY_NO;
  
  /** was polarity switching used*/
  private boolean polaritySwitching_ = false;
  
  /**
   * Constructs an AbstractXmlSpectraReader object with given information.
   * 
   * @param callbacks
   *          Pass an array of objects implementing this interface. The methods will be
   *          called whenever a msScan header or individual scans are generated
   *          during the reading process.
   * @param parseMsMS shall the MS/MS spectra be parsed too
   */
  public AbstractXMLSpectraReader(AddScan[] callbacks, boolean parseMsMs)
  {
    this.adders_ = callbacks;
    this.generateMzThresholds();
    this.parseMsMs_ = parseMsMs;
  }
  
  /**
   * Constructs an AbstractXmlSpectraReader object with given information.
   * 
   * @param callbacks
   *          Pass an array of objects implementing this interface. The methods will be
   *          called whenever a msScan header or individual scans are generated
   *          during the reading process.
   * @param parseMsMs shall the MS/MS spectra be parsed too
   * @param multiplicationFactorForInt the multiplication factor to be used, to create integer values out of the m/z float values
   */
  public AbstractXMLSpectraReader(AddScan[] callbacks, boolean parseMsMs, int multiplicationFactorForInt)
  {
    this(callbacks, parseMsMs);
    this.multiplicationFactorForInt_ = multiplicationFactorForInt;
    this.lowestMz_ = ONE_MILLION * this.multiplicationFactorForInt_;  
  }
  
    
  /**
   * Method to read in an XML file synchronously. This method returns after
   * reading in the complete XML file.
   * 
   * @param fileName
   *          File Name (including full path if required) for the XML file to
   *          be read
   * @throws CgException
   *           All internal exceptions are mapped to the CgException type.
   */
  public void ReadFile(String fileName) throws CgException
  {
    this.ReadFile(fileName,false); 
  }

  /**
   * Method to read in an XML file synchronously. This method returns after
   * reading in the complete XML or a defined subset. If the required XML start
   * tag is not found, it throws a CgException.
   * 
   * @param fileName
   *          File Name (including full path if required) for the XML file to
   *          be read
   * @param readOnlyRequiredInfoForMultiThreading 
   *          Reads only the information required for getting an idea about the 
   *          m/z ranges for splitting in iterations and/or threads,
   *          the retention time range, 
   *          highest msLevel, 
   *          as well as whether polarity switching is used
   * @throws CgException
   *          All internal exceptions are mapped to the CgException type.
   */
  public void ReadFile(String fileName, boolean readOnlyRequiredInfoForMultiThreading) throws CgException
  {
    int eventType;
    boolean foundTagRun = false;
    
    try {
      
      // =========================================================
      // Open the file:
      // =========================================================

      XMLInputFactory factory = XMLInputFactory.newInstance();
      factory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, true);
      factory.setProperty(XMLInputFactory.IS_COALESCING, true);
      File XMLFile = new File(fileName);

      if (XMLFile.exists()){
        inStream_ = new FileInputStream(XMLFile); 
      }else{
        XMLFile = new File(fileName+".gz");
        if (XMLFile.exists()){
          inStream_ = new GZIPInputStream(new FileInputStream(XMLFile)); 
        }else{
          throw new CgException(String.format("The file %s does not exist!", fileName));
        }
      }
      reader_ = factory.createXMLStreamReader(inStream_);
      
      // =========================================================
      // Read the XML Data:
      // =========================================================
      
      eventType = reader_.getEventType();
      while (reader_.hasNext()) {
        switch (eventType) {
          case XMLStreamReader.START_DOCUMENT:
            break;
          case XMLStreamReader.END_DOCUMENT:
            break;
          case XMLStreamReader.START_ELEMENT:
            if (reader_.getLocalName().equalsIgnoreCase(getTagRun())) {
              foundTagRun = true;
              readMsRun(readOnlyRequiredInfoForMultiThreading);
            }
            break;
          default:
            break;
        }
        eventType = reader_.next();
      }
    }
    catch (Exception ex) {
      ex.printStackTrace();
      throw new CgException(ex.getMessage());
    }finally{
      if (!foundTagRun) {
        throw new CgException(String.format(
            "The file %s does not contain the required tag <%s> and could therefore not be read by %s.", 
            fileName, getTagRun(), this.getClass().getName()));
      }
      try{
        inStream_.close();
        reader_ = null;
      }catch(IOException iox){
        iox.printStackTrace();
      }
    }
  }
  
  /**
   * This method gets the value of the required attribute from reader and throws an exception if 
   * the attribute is not found.
   * 
   * @param reader XMLStreamReader instance used for parsing
   * @param attribute Attribute to be found
   * @return the value of the attribute
   */
  public String getRequiredAttribute(XMLStreamReader reader, String attribute) {
    String attributeValue = reader.getAttributeValue(null, attribute);
    if (attributeValue == null) {
      throw new IllegalStateException(String.format("Tag %s must provide the attribute `%s`(Line %s)", 
          reader.getLocalName(), attribute, reader.getLocation().getLineNumber()));
    }
    return attributeValue;
  }
  
  /**
   * If the given values are currently maxima in the file,
   * the global values are set to the given values.
   * 
   * @param lowMz The lowest m/z value found.
   * @param highMz The highest m/z value found.
   */
  protected void setCurrentGlobalMaxima(float lowMz, float highMz) {
    int currentLowMz =  Math.round(lowMz*getMultiplicationFactorForInt());
    int currentHighMz = Math.round(highMz*getMultiplicationFactorForInt());                
    if (currentLowMz<getLowestMz()) setLowestMz(currentLowMz);
    if (currentHighMz>getHighestMz()) setHighestMz(currentHighMz);
  }
  
  /**
   * This method uncompresses a zlib compressed byte array
   * @param decoded The compressed byte array
   * 
   * @return the uncompressed byte array
   */
  protected static byte[] decompressZLIB(byte[] decoded) {
    Inflater decompressor = new Inflater();
    decompressor.setInput(decoded);
    try (ByteArrayOutputStream bos = new ByteArrayOutputStream(decoded.length)) {
        byte[] buf = new byte[1024];
        while (!decompressor.finished()) {
          int count = decompressor.inflate(buf);
          bos.write(buf, 0, count);
        }
        decoded = bos.toByteArray();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
    return decoded;
  }
  
  /**
   * This method generates the maxRange object from the individual AddScan ranges.
   */
  private void generateMzThresholds(){
    float lowest = Float.MAX_VALUE;
    float highest = 0;
    for (int i=0; i!=this.adders_.length; i++){
      AddScan adder = this.adders_[i];
      if (adder.getLowerThreshold()<lowest) lowest = adder.getLowerThreshold();
      if (adder.getUpperThreshold()>highest) highest = adder.getUpperThreshold();
    }
    this.maxRange_ = new Range(lowest,highest);
  }
  
  /**
   * This method generates an "empty space" separated string containing the precursor m/z values of the current spectrum
   * 
   * @param precursorMzs The single precursor m/z values
   * 
   * @return the "empty space" separated string containing the precursor m/z values of the current spectrum
   */
  protected String getPrecursorMzString (Vector<String> precursorMzs){
    if (precursorMzs.size()==1) return precursorMzs.get(0);
    StringBuilder bd = new StringBuilder();
    for (int i=0 ; i!=precursorMzs.size(); i++){
      bd.append(precursorMzs.get(i));
      if (i!=(precursorMzs.size()-1)) bd.append(" ");
    }
    return bd.toString();
  }
  
  /**
   * This method reads the full structure of the XML element constant TAG_RUN into memory. 
   * Basically it reads each scan with msLevel=1, inner scans are registered, but not fully read.
   * 
   * @param readOnlyRequiredInfoForMultiThreading
   *          Reads only the information required for getting an idea about the 
   *          m/z ranges for splitting in iterations and/or threads,
   *          the retention time range, 
   *          highest msLevel, 
   *          as well as whether polarity switching is used
   * 
   * @throws CgException 
   *          All internal exceptions are mapped to the CgException type.
   */
  protected abstract void readMsRun(boolean readOnlyRequiredInfoForMultiThreading) throws CgException;
  
  /**
   * This method reads only the information required for getting an idea about the 
   * m/z ranges for splitting in iterations and/or threads,
   * the retention time range, 
   * highest msLevel, 
   * as well as whether polarity switching is used
   *          
   * @throws CgException All internal exceptions are mapped to the CgException type.
   */
  protected abstract void readOnlyRequiredInfoForMultiThreading() throws CgException;
  
  /**
   * Returns the XML element constant, the content of which will be read. 
   * Everything outside of it will be ignored by the reader.
   * 
   * @return the XML element constant
   */
  public abstract String getTagRun();
  
  /**
   * Returns the maximally allowed m/z range provided by all the AddScan interfaces
   * 
   * @return the maximally allowed m/z range
   */
  protected Range getMaxRange()
  {
    return this.maxRange_;
  }
  
  /**
   * Sets the AddScan interfaces for this reader
   * 
   * @param adders The AddScan interfaces to be used by this reader
   */
  public void setAdders(AddScan[] adders)
  {
    this.adders_ = adders;
  }
  
  /**
   * Returns a boolean which defines whether MS/MS spectra are parsed too
   * 
   * @return boolean which defines whether MS/MS spectra are parsed too
   */
  protected boolean getParseMsMs()
  {
    return this.parseMsMs_;
  }
  
  /**
   * Returns the multiplication factor to be used, to create integer values out of the m/z float values
   * 
   * @return an integer representing the multiplication factor
   */
  protected int getMultiplicationFactorForInt()
  {
    return this.multiplicationFactorForInt_;
  }
  
  /**
   * Sets the highest available m/z value
   * 
   * @param highestMz The highest available m/z value
   */
  protected void setHighestMz(int highestMz)
  {
    this.highestMz_ = highestMz;
  }
  
  /**
   * Returns the highest available m/z value
   * 
   * @return the highest available m/z value
   */
  public int getHighestMz()
  {
    return this.highestMz_;
  }
  
  /**
   * Sets the lowest available m/z value
   * 
   * @param lowestMz the lowest available m/z value
   */
  protected void setLowestMz(int lowestMz)
  {
    this.lowestMz_ = lowestMz;
  }

  /**
   * Returns the lowest available m/z value
   * 
   * @return the lowest available m/z value
   */
  public int getLowestMz()
  {
    return this.lowestMz_;
  }
  
  /**
   * Sets a boolean which defines whether polarity switching is used
   * 
   * @param polaritySwitching whether polarity switching is used
   */
  protected void setPolaritySwitching(boolean polaritySwitching)
  {
    this.polaritySwitching_ = polaritySwitching;
  }
  
  /**
   * Returns a boolean which defines whether polarity switching is used
   * 
   * @return boolean which defines whether polarity switching is used
   */
  public boolean usesPolaritySwitching()
  {
    return this.polaritySwitching_;
  }
  
  /** @deprecated*/
  //this method is obsolete; only there because of the interface
  public void setLowerThreshold(float lowerThreshold)
  {
  }

  /** @deprecated*/
  //this method is obsolete; only there because of the interface
  public void setUpperThreshold(float upperThreshold)
  {
  }
  
}
