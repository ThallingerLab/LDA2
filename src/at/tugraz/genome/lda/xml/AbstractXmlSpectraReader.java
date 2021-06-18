/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

//import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
//import java.nio.ByteBuffer;
//import java.nio.DoubleBuffer;
//import java.nio.FloatBuffer;
//import java.util.Hashtable;
import java.util.Vector;
import java.util.zip.GZIPInputStream;
//import java.util.zip.Inflater;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;

//TODO: the following lines are for QQQ PIS/NLS - activate when necessary
//import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.swing.Range;
//import at.tugraz.genome.maspectras.quantification.CgBase64;
import at.tugraz.genome.maspectras.quantification.CgDefines;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgScan;
import at.tugraz.genome.maspectras.quantification.CgScanHeader;
//import at.tugraz.genome.maspectras.quantification.MsMsScan;
//import at.tugraz.genome.maspectras.utils.StringUtils;

/**
 * StAX based mzXML Reader for the Package.
 * 
 * This class reads an mzXML file with hierarchical scans. It uses two methods
 * defined in the AddScan interface to return created header and spectrum
 * information to the caller.
 * 
 * @author Juergen Hartler
 * @author Leonida Lamp
 */
public abstract class AbstractXmlSpectraReader implements XmlSpectraReader
{
  //the reader
  XMLStreamReader rdr_;
  
  //the reader can serve several AddScan interfaces
  AddScan[] adders_;
  
  // if true MS/MS spectra are parsed too
  protected boolean parseMsMs_;
  
  // the lowest available m/z value in integer format (original m/z value times the multiplication factor for int)
  protected int lowestMz_ = 1000000 * CgDefines.mzMultiplicationFactorForInt;
  // the highest available m/z value in integer format (original m/z value times the multiplication factor for int)
  protected int highestMz_ = 0;
  
  // the provided input stream
  private InputStream inStream_;
  // the multiplication factor to be used, to create integer values out of the m/z float values
  protected int multiplicationFactorForInt_;
  
  //the call back interface for the header information
  protected CgScanHeader myHeader_;
  
  //the maximally allowed m/z range provided by all the AddScan interfaces
  protected Range maxRange_ = null;
  
  /** the precursor m/z values of the current spectrum*/
  protected Vector<String> precursorMz_;
  
  /** the MS-level of the previous scan*/
  protected int lastMsLevel_;
  
  /** the polarity of the current scan*/
  protected int currentPolarity_ = CgDefines.POLARITY_NO;
  
  /** was polarity switching used*/
  protected boolean polaritySwitching_ = false;
  
  /**
   * @param callbacks
   *          Pass an array of objects implementing this interface. The methods will be
   *          called whenever a msScan header or individual scans are generated
   *          during the reading process.
   * @param parseMsMS shall the MS/MS spectra be parsed too
   */
  public AbstractXmlSpectraReader(AddScan[] callbacks, boolean parseMsMs)
  {
    multiplicationFactorForInt_ = CgDefines.mzMultiplicationFactorForInt;
    adders_ = callbacks;
    this.parseMsMs_ = parseMsMs;
    lowestMz_ = 1000000 * CgDefines.mzMultiplicationFactorForInt;
    currentPolarity_ = CgDefines.POLARITY_NO;
    polaritySwitching_ = false;
  }
  
  /**
   * 
   * @param callbacks
   *          Pass an array of objects implementing this interface. The methods will be
   *          called whenever a msScan header or individual scans are generated
   *          during the reading process.
   * @param parseMsMs shall the MS/MS spectra be parsed too
   * @param multiplicationFactorForInt the multiplication factor to be used, to create integer values out of the m/z float values
   */
  public AbstractXmlSpectraReader(AddScan[] callbacks, boolean parseMsMs, int multiplicationFactorForInt)
  {
    this(callbacks, parseMsMs);
    multiplicationFactorForInt_ = multiplicationFactorForInt;
    lowestMz_ = 1000000 * multiplicationFactorForInt_;  
  }
  
    
  /**
   * Method to read in an mzXML file synchronously. This method returns after
   * reading in the complete mzXML file.
   * 
   * @param fileName
   *          File Name (including full path if required) for the mzXML file to
   *          be read
   * @throws CgException
   *           All internal exceptions are mapped to the CgException type.
   */
  public void ReadFile(String fileName) throws CgException
  {
    this.ReadFile(fileName,false); 
  }

  /**
   * Method to read in an mzXML file synchronously. This method returns after
   * reading in the complete mzXML or a defined subset.
   * 
   * @param fileName
   *          File Name (including full path if required) for the mzXML file to
   *          be read
   * @param readJustMzMaxima reads only the m/z maxima - for getting an idea
   *          about the m/z ranges for splitting in iterations and/or threads
   * @throws CgException
   *           All internal exceptions are mapped to the CgException type.
   */
  public void ReadFile(String fileName, boolean readJustMzMaxima) throws CgException
  {
    this.generateMzThresholds();
    
    int eventType;
    try {
      // =========================================================
      // Open the file:
      // =========================================================

      XMLInputFactory factory = XMLInputFactory.newInstance();
      factory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, true);
      factory.setProperty(XMLInputFactory.IS_COALESCING, true);
      File mzXMLFile = new File(fileName);

      if (mzXMLFile.exists()){
        inStream_ = new FileInputStream(mzXMLFile); 
      }else{
        mzXMLFile = new File(fileName+".gz");
        if (mzXMLFile.exists()){
          inStream_ = new GZIPInputStream(new FileInputStream(mzXMLFile)); 
        }else{
          throw new CgException("The file "+fileName+" does not exist!");
        }
      }
      rdr_ = factory.createXMLStreamReader(inStream_);
      // =========================================================
      // Read the XML Data:
      // =========================================================
      eventType = rdr_.getEventType();
      do {
        switch (eventType) {
          case XMLStreamReader.START_DOCUMENT:
            break;
          case XMLStreamReader.END_DOCUMENT:
            break;
          case XMLStreamReader.START_ELEMENT:
            /**
             * TODO: attempt to make this if statement more general, in mzXML: "msRun", in mzML: "run", can probably be generalized to "run"
             */
            if (rdr_.getLocalName().equalsIgnoreCase("msRun"))
              XmlReadMsRun(readJustMzMaxima);
            break;
        }
        eventType = rdr_.next();
      } while (eventType != XMLStreamReader.END_DOCUMENT);
    }
    catch (Exception ex) {
      ex.printStackTrace();
      throw new CgException(ex.getMessage());
    }finally{
      try{
        inStream_.close();
        rdr_ = null;
      }catch(IOException iox){
        iox.printStackTrace();
      }
    }
  }
  
  /**
   * This Function reads a full MsRun Structure into Memory. Basically it
   * reads each scan with msLevel=1, inner scans are registered, but not fully
   * read.
   * 
   * @param readJustMzMaxima reads only the m/z maxima - for getting an idea
   *          about the m/z ranges for splitting in iterations and/or threads
   * 
   * @throws CgException
   */
  protected abstract void XmlReadMsRun(boolean readJustMzMaxima) throws CgException;
  
  /**
   * this function reads only the m/z maxima - for getting an idea
   *      about the m/z ranges for splitting in iterations and/or threads
   *          
   * @throws CgException
   */
  protected abstract void xmlReadMaxima() throws CgException;
  
  /**
   * this function reads out the m/z maxima of a <peaks> tag - for getting an idea
   *      about the m/z ranges for splitting in iterations and/or threads
   *      
   * @param peaksCount the number of peaks in the <peaks> tag
   * @return
   * @throws CgException
   */
  protected abstract float[] readMaximaFromPeaks(int peaksCount) throws CgException;
  
  /**
   * This Method reads a scan. It calls itself recursively in case scans are
   * element of our scan. However, if the scan level is > 1, the scan's content
   * is skipped for performance reasons.
   * 
   * @param scBase1
   *          Pass the CgScan object that represents the level 1 scan to which
   *          this scan belongs to.
   * @param ranges1 the individually allowed m/z ranges of the AddScan interfaces
   * @param maxRange the maximally allowed m/z range provided by all the AddScan interface
   * 
   * @throws CgException
   */
  protected abstract void XmlReadScan(Vector<CgScan> scBase1, Vector<Range> ranges1, Range maxRange) throws CgException;
  

  /**
   * This Function reads the peaks of a single scan. If there is a valid CgScan
   * passed as parameter, it stores the values. If not, values are just skipped.
   * 
   * @param scans
   *          Pass the CgScan object that represents the level 1 scan to which
   *          this scan belongs to.
   * @param ranges
   *          The m/z range restrictions that apply for every CgScan object
   * @param maxRange the maximally allowed m/z range provided by all the AddScan interface
   * @param peaksCount the number of peaks in the <peaks> tag
   * @param msms true if it is an MSn scan
   * @param foundMzBorders true if it has to be scanned for the lowest and highest m/z values of the file
   * 
   * @throws CgException
   */
  protected abstract void XmlReadPeaks(Vector<CgScan> scans, Vector<Range> ranges, Range maxRange, int peaksCount, boolean msms, boolean foundMzBorders) throws CgException;
  
  /**
   * Function converts a Time Interval to a double value[s]. Do something
   * better, this function is really very basic.
   * 
   * @param s
   *          String to parse
   * @return float representing a number of seconds
   */
  public static float XmlTimeIntervalToTime(String s)
  {
    if (s.startsWith("PT"))
      s = s.substring(2);
    if (s.endsWith("S"))
      s = s.substring(0, s.length() - 1);
    s = s.replace('e', 'E');
    return Float.parseFloat(s);
  }

  public int getHighestMz()
  {
    return this.highestMz_;
  }

  public int getLowestMz()
  {
    return this.lowestMz_;
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
  
  /**
   * sets the AddScan interfaces for this reader
   * @param adders the AddScan interfaces to be used by this reader
   */
  public void setAdders(AddScan[] adders)
  {
    this.adders_ = adders;
  }
  
  /**
   * generates the maxRange object from the individual AddScan ranges
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
   * generates an "empty space" separated string containing the precursor m/z values of the current spectrum
   * @param precursorMzs the single precursor m/z values
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

  public boolean usesPolaritySwitching()
  {
    return this.polaritySwitching_;
  }

}
