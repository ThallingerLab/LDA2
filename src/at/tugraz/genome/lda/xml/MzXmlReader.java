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

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.util.Hashtable;
import java.util.Vector;
import java.util.zip.GZIPInputStream;
import java.util.zip.Inflater;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;

//TODO: the following lines are for QQQ PIS/NLS - activate when necessary
/****import at.tugraz.genome.lda.LipidomicsConstants;*/
import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.maspectras.quantification.CgBase64;
import at.tugraz.genome.maspectras.quantification.CgDefines;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgScan;
import at.tugraz.genome.maspectras.quantification.CgScanHeader;
import at.tugraz.genome.maspectras.quantification.MsMsScan;
import at.tugraz.genome.maspectras.utils.StringUtils;

/**
 * StAX based mzXML Reader for the Package.
 * 
 * This class reads an mzXML file with hierarchical scans. It uses two methods
 * defined in the AddScan interface to return created header and spectrum
 * information to the caller.
 * 
 * @author Juergen Hartler
 */
public class MzXmlReader implements XmlSpectraReader
{
  //the reader
  XMLStreamReader rdr_;
  
  //the reader can serve several AddScan interfaces
  AddScan[] adders_;
  
  // if true MS/MS spectra are parsed too
  private boolean parseMsMs_;
  
  // the lowest available m/z value in integer format (original m/z value times the multiplication factor for int)
  private int lowestMz_ = 1000000 * CgDefines.mzMultiplicationFactorForInt;
  // the highest available m/z value in integer format (original m/z value times the multiplication factor for int)
  private int highestMz_ = 0;
  
  // the provided input stream
  private InputStream inStream_;
  // the multiplication factor to be used, to create integer values out of the m/z float values
  private int multiplicationFactorForInt_;
  
  //the call back interface for the header information
  private CgScanHeader myHeader_;
  
  //the maximally allowed m/z range provided by all the AddScan interfaces
  private Range maxRange_ = null;
  
  /** the precursor m/z values of the current spectrum*/
  private Vector<String> precursorMz_;
  
  /** the MS-level of the previous scan*/
  private int lastMsLevel_;
  
  /** the polarity of the current scan*/
  private int currentPolarity_ = CgDefines.POLARITY_NO;
  
  /** was polarity switching used*/
  private boolean polaritySwitching_ = false;
  
  /**
   * @param callbacks
   *          Pass an array of objects implementing this interface. The methods will be
   *          called whenever a msScan header or individual scans are generated
   *          during the reading process.
   * @param parseMsMS shall the MS/MS spectra be parsed too
   */
  public MzXmlReader(AddScan[] callbacks, boolean parseMsMs)
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
  public MzXmlReader(AddScan[] callbacks, boolean parseMsMs, int multiplicationFactorForInt)
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
   * reads each scan with msLevel=1, inner scans are re- gistered, but not fully
   * read.
   * 
   * @param readJustMzMaxima reads only the m/z maxima - for getting an idea
   *          about the m/z ranges for splitting in iterations and/or threads
   * 
   * @throws CgException
   */
  private void XmlReadMsRun(boolean readJustMzMaxima) throws CgException
  {
    int i;
    int eventType;
    lastMsLevel_ = 0;
    
    myHeader_ = new CgScanHeader();
    myHeader_.highestMSLevel = 1;
    // =========================================================
    // Read the msRun-Attributes. We put these into our Header
    // object.
    // =========================================================

    for (i = 0; i < rdr_.getAttributeCount(); i++) {
      if (rdr_.getAttributeLocalName(i).equalsIgnoreCase("scanCount")) {
        myHeader_.ScanCount = Integer.parseInt(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "startTime") {
        myHeader_.StartTime = XmlTimeIntervalToTime(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "endTime") {
        myHeader_.EndTime = XmlTimeIntervalToTime(rdr_.getAttributeValue(i));
      }
    }
    if (adders_ != null && adders_.length>0)
      for (AddScan adder : adders_) adder.AddHeader(myHeader_);
    else{
      if (!readJustMzMaxima)
        throw new CgException("No adder for Header and Scans defined.");
    }  

    // =========================================================
    // Finally read the bottom level scans.
    // =========================================================
    try {
      do {
        eventType = rdr_.getEventType();
        switch (eventType) {
          case XMLStreamReader.START_ELEMENT:
            if (rdr_.getLocalName().equalsIgnoreCase("parentFile")){
              for (i = 0; i < rdr_.getAttributeCount(); i++) {
                if (rdr_.getAttributeLocalName(i).equalsIgnoreCase("fileName")) {
                  for (AddScan adder : adders_) adder.addParentFileName(StringUtils.getFileNameWOSuffix(rdr_.getAttributeValue(i)));
                }
              }
            }else if (rdr_.getLocalName().equalsIgnoreCase("scan")){
              if (readJustMzMaxima)
                xmlReadMaxima();
              else  
                XmlReadScan(null,null,this.maxRange_);
            }  
            break;
          case XMLStreamReader.END_ELEMENT:
            if (rdr_.getLocalName().equalsIgnoreCase("msRun"))
              return;
            break;
        }
        eventType = rdr_.next();
      } while (eventType != XMLStreamReader.END_DOCUMENT);
    }
    catch (Exception ex) {
      ex.printStackTrace();
      throw new CgException(ex.getMessage());
    }
  }
  
  /**
   * this function reads only the m/z maxima - for getting an idea
   *      about the m/z ranges for splitting in iterations and/or threads
   *          
   * @throws CgException
   */
  private void xmlReadMaxima() throws CgException
  {
    int i;
    int eventType;

    int num = -1;
    int msLevel = 0;
    float lowMz = 0;
    float highMz = 0;
    boolean lowMzFound = false;
    boolean highMzFound = false;
    String polarityString = "";
    int polarity = CgDefines.POLARITY_NO;
    int peaksCount = 0;
    //TODO: the following lines are for QQQ PIS/NLS - activate when necessary
/***    float precursor = -1;*/

    for (i = 0; i < rdr_.getAttributeCount(); i++) {
      if (rdr_.getAttributeLocalName(i) == "num") {
        num = Integer.parseInt(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "msLevel") {
        msLevel = Integer.parseInt(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "lowMz") {
        lowMz = Float.parseFloat(rdr_.getAttributeValue(i));
        lowMzFound = true;
      } else if (rdr_.getAttributeLocalName(i) == "highMz") {
        highMz = Float.parseFloat(rdr_.getAttributeValue(i));
        highMzFound=true;
      } else if (rdr_.getAttributeLocalName(i) == "peaksCount") {
        peaksCount = Integer.parseInt(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "polarity") {
        polarityString = rdr_.getAttributeValue(i);
        if (polarityString.equalsIgnoreCase("+"))
          polarity = CgDefines.POLARITY_POSITIVE;
        else if (polarityString.equalsIgnoreCase("-"))
          polarity = CgDefines.POLARITY_NEGATIVE;
        else
          throw new CgException("The scan contains an unknown polarity \""+polarityString+"\" at scan number: "+num);
      }
    }
    boolean foundMzBorders = false;
    if (lowMzFound&&highMzFound) foundMzBorders = true;
    else if (msLevel==1 && !myHeader_.hasMS1Scans){
      lowestMz_ = 1000000 * CgDefines.mzMultiplicationFactorForInt;
      highestMz_ = 0;      
    }
    if (msLevel==1)
      myHeader_.hasMS1Scans = true;

    if (polarity!=CgDefines.POLARITY_NO){
      if (currentPolarity_==CgDefines.POLARITY_NO)
        currentPolarity_ = polarity;
      // if the polarity is suddenly different, polarity switching is used
      else if (currentPolarity_!=polarity)
        polaritySwitching_ = true;
    }
    
    try {
      eventType = rdr_.next();
      do {
        switch (eventType) {
          case XMLStreamReader.START_ELEMENT:
            if (rdr_.getLocalName() == "peaks") {
              if (msLevel == 1) {
                if (!foundMzBorders){
                  float[] maxima = readMaximaFromPeaks(peaksCount);
                  lowMz = maxima[0];
                  highMz = maxima[1];
                }
                int currentLowMz =  Math.round(lowMz*multiplicationFactorForInt_);
                int currentHighMz = Math.round(highMz*multiplicationFactorForInt_);                
                if (currentLowMz<this.lowestMz_) this.lowestMz_ = currentLowMz ;
                if (currentHighMz>this.highestMz_) this.highestMz_ = currentHighMz;
                // =================================================
                // Now we can inform our caller that we have a valid
                // new scan!
                // =================================================
                
                //TODO: the following lines are for QQQ PIS/NLS - activate when necessary
/****              } else if (!myHeader_.hasMS1Scans && LipidomicsConstants.isShotgun() && precursor>=0){
                int currentLowMz =  Math.round((precursor*0.999f)*multiplicationFactorForInt_);
                int currentHighMz = Math.round((precursor*1.001f)*multiplicationFactorForInt_);                
                if (currentLowMz<this.lowestMz_) this.lowestMz_ = currentLowMz ;
                if (currentHighMz>this.highestMz_) this.highestMz_ = currentHighMz;*/
              }
              if (parseMsMs_ && msLevel > myHeader_.highestMSLevel) myHeader_.highestMSLevel=msLevel;
            //TODO: the following lines are for QQQ PIS/NLS - activate when necessary
/****            } else if (rdr_.getLocalName().equalsIgnoreCase("precursorMz") && !myHeader_.hasMS1Scans && LipidomicsConstants.isShotgun()) {
              int attributeCount = rdr_.getAttributeCount();
              for (i = 0; i < attributeCount; i++) {
                if (rdr_.getAttributeLocalName(i) == "precursorIntensity") {
                  try {
                    rdr_.next();
                    String childNode = rdr_.getText().trim();
                    if (childNode != null) {
                      precursor = Float.parseFloat(childNode);
                      break;
                    }
                  }
                  catch (Exception ex) {
                    throw new CgException(ex.getMessage());
                  }
                }
              }*/
            } else if (rdr_.getLocalName().equalsIgnoreCase("scan")) {
              xmlReadMaxima();
            }
            break;
          case XMLStreamReader.END_ELEMENT:
            if (rdr_.getLocalName().equalsIgnoreCase("scan"))
              return;
            break;
        }
        eventType = rdr_.next();
      } while (eventType != XMLStreamReader.END_DOCUMENT);
    }
    catch (Exception ex) {
      throw new CgException(ex.getMessage());
    }
    
  }
  
  /**
   * this function reads out the m/z maxima of a <peaks> tag - for getting an idea
   *      about the m/z ranges for splitting in iterations and/or threads
   *      
   * @param peaksCount the number of peaks in the <peaks> tag
   * @return
   * @throws CgException
   */
  private float[] readMaximaFromPeaks(int peaksCount) throws CgException{
    float[] maxima = new float[]{Float.MAX_VALUE,0f};
    if (peaksCount<1) return maxima;
    int precision = -1;
    String compressionType = null;
    CgBase64 cgb = new CgBase64();
    for (int i = 0; i < rdr_.getAttributeCount(); i++) {
      if (rdr_.getAttributeLocalName(i).equalsIgnoreCase("precision")) {
        precision = Integer.parseInt(rdr_.getAttributeValue(i));
      }  else if (rdr_.getAttributeLocalName(i) == "compressionType") {
        compressionType = rdr_.getAttributeLocalName(i);
      }
    }
    try {
      rdr_.next();
      String s = rdr_.getText().trim(); // In s we have a Base64 coded value array!
      if (s == null || s.equalsIgnoreCase("</peaks>"))
        return maxima;

      // =================================================
      // Process the data, in case we have to store it. In
      // C#, we have to store the byte array into a memory
      // stream and later on read it out float by float by
      // a binary reader.
      // =================================================
      byte[] decoded = cgb.decode(s);
      if (compressionType!=null && compressionType.equalsIgnoreCase("zlib")){
        Inflater decompressor = new Inflater();
        decompressor.setInput(decoded);
        ByteArrayOutputStream bos = null;
        try {
            bos = new ByteArrayOutputStream(decoded.length);
            byte[] buf = new byte[1024];
            while (!decompressor.finished()) {
                int count = decompressor.inflate(buf);
                bos.write(buf, 0, count);
            }

        } finally {
            try {
                bos.close();
            } catch (Exception nope) { /* This exception doesn't matter */ }
        }
        decoded = bos.toByteArray();
      }
      
      ByteBuffer byteBuf = ByteBuffer.wrap(decoded);
      float lowestMzValue = Float.MAX_VALUE;
      float highestMzValue = 0f;
      if (precision==64){
        double doubleArray[] = new double[decoded.length/8];
        DoubleBuffer doubleBuf = byteBuf.asDoubleBuffer();
        doubleBuf.get(doubleArray);
        lowestMzValue = (float)doubleArray[0];
        highestMzValue = (float)doubleArray[(peaksCount-1)*2];
      }else{
        float floatArray[] = new float[decoded.length / 4];
        FloatBuffer floatBuf = byteBuf.asFloatBuffer();
        floatBuf.get(floatArray);
        lowestMzValue = floatArray[0];
        highestMzValue = floatArray[(peaksCount-1)*2];
      }
      maxima[0] = lowestMzValue;
      maxima[1] = highestMzValue;
    }
    catch (Exception ex) {
      ex.printStackTrace();
      throw new CgException(ex.getMessage());
    }
    return maxima;
  }
  
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
  private void XmlReadScan(Vector<CgScan> scBase1, Vector<Range> ranges1, Range maxRange) throws CgException
  {
    int i;
    int eventType;
    Vector<CgScan> baseScans = scBase1;
    Vector<Range> scanRanges = ranges1;
    CgScan sc = null;
    int num = 0;
    int msLevel = 0;
    int peaksCount = 0;
    float retentionTime = 0;
    float lowMz = 0;
    float highMz = 0;
    float basePeakMz = 0;
    float basePeakIntensity = 0;
    float totIonCurrent = 0;
    float precursorIntensity = 0;
    String polarityString = "";
    int polarity = CgDefines.POLARITY_NO;
    boolean lowMzFound = false;
    boolean highMzFound = false;
    //TODO: the following lines are for QQQ PIS/NLS - activate when necessary
    /****float precursorMz = -1f;*/

    // =========================================================
    // First of all we read the attributes:
    // =========================================================

    for (i = 0; i < rdr_.getAttributeCount(); i++) {
      if (rdr_.getAttributeLocalName(i) == "num") {
        num = Integer.parseInt(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "msLevel") {
        msLevel = Integer.parseInt(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "peaksCount") {
        peaksCount = Integer.parseInt(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "retentionTime") {
        retentionTime = XmlTimeIntervalToTime(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "lowMz") {
        lowMz = Float.parseFloat(rdr_.getAttributeValue(i));
        lowMzFound = true;
      } else if (rdr_.getAttributeLocalName(i) == "highMz") {
        highMz = Float.parseFloat(rdr_.getAttributeValue(i));
        highMzFound = true;
      } else if (rdr_.getAttributeLocalName(i) == "basePeakMz") {
        basePeakMz = Float.parseFloat(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "basePeakIntensity") {
        basePeakIntensity = Float.parseFloat(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "totIonCurrent") {
        totIonCurrent = Float.parseFloat(rdr_.getAttributeValue(i));
      } else if (rdr_.getAttributeLocalName(i) == "polarity") {
        polarityString = rdr_.getAttributeValue(i);
        if (polarityString.equalsIgnoreCase("+"))
          polarity = CgDefines.POLARITY_POSITIVE;
        else if (polarityString.equalsIgnoreCase("-"))
          polarity = CgDefines.POLARITY_NEGATIVE;
        else
          throw new CgException("The scan contains an unknown polarity \""+polarityString+"\" at scan number: "+num);
      }
    }
    boolean foundMzBorders = false;
    if (lowMzFound&&highMzFound) foundMzBorders = true;
    else if (msLevel==1 && !myHeader_.hasMS1Scans){
      lowestMz_ = 1000000 * CgDefines.mzMultiplicationFactorForInt;
      highestMz_ = 0;      
    }
    if (msLevel==1)
      myHeader_.hasMS1Scans = true;
    if (polarity!=CgDefines.POLARITY_NO){
      if (currentPolarity_==CgDefines.POLARITY_NO)
        currentPolarity_ = polarity;
      // if the polarity is suddenly different, polarity switching is used
      else if (currentPolarity_!=polarity)
        polaritySwitching_ = true;
    }

    
    if (msLevel<3){
      precursorMz_ = new Vector<String>();
    } else if (lastMsLevel_>=msLevel){
      for (i=(lastMsLevel_+1); i!=msLevel; i--){
        precursorMz_.remove(precursorMz_.size()-1);
      }
    }
    lastMsLevel_ = msLevel;
    
    // =========================================================
    // Now we read the peaks:
    // =========================================================

    try {
      eventType = rdr_.next();
      do {
        switch (eventType) {
          case XMLStreamReader.START_ELEMENT:
            if (rdr_.getLocalName().equalsIgnoreCase("peaks")) {
              if (msLevel == 1) {
                if (msLevel > myHeader_.highestMSLevel) myHeader_.highestMSLevel=msLevel;
                //sc = new CgScan(peaksCount);
                sc = new CgScan(0);
                sc.Num = num;
                sc.MsLevel = msLevel;
                sc.RetentionTime = retentionTime;
                sc.LowMz = lowMz;
                sc.HighMz = highMz;
                sc.BasePeakMz = basePeakMz;
                sc.BasePeakIntensity = basePeakIntensity;
                sc.TotIonCurrent = totIonCurrent;
                sc.setPolarity(polarity);
                if (adders_ != null && adders_.length>0){
                  Vector<CgScan> scans = new Vector<CgScan>();
                  Vector<Range> ranges = new Vector<Range>();
                  for (AddScan adder : adders_){
                    if (adders_.length==1) scans.add(sc);
                    else scans.add(new CgScan(sc));
                    ranges.add(new Range(adder.getLowerThreshold(),adder.getUpperThreshold()));
                  }

                  XmlReadPeaks(scans,ranges,maxRange,peaksCount,false,foundMzBorders);
                  for (int j=0; j!=adders_.length; j++){
                    AddScan adder = adders_[j];
                    sc = scans.get(j);
                    int currentLowMz = Math.round(sc.LowMz*multiplicationFactorForInt_);
                    int currentHighMz = Math.round(sc.HighMz*multiplicationFactorForInt_);
                    if (currentLowMz<this.lowestMz_) this.lowestMz_ = currentLowMz ;
                    if (currentHighMz>this.highestMz_) this.highestMz_ = currentHighMz;
                    // =================================================
                    // Now we can inform our caller that we have a valid
                    // new scan!
                    // =================================================
                    adder.AddScan(sc);
                  }
                  if (msLevel > myHeader_.highestMSLevel) myHeader_.highestMSLevel=msLevel;
                } else
                  throw new CgException(
                      "No adder for Header and Scans defined.");
              } else {
                //TODO: the following lines are for QQQ PIS/NLS - activate when necessary
                //this generates an artificial MS1 scans for shotgun MSn-only data
/****                if (!myHeader_.hasMS1Scans && LipidomicsConstants.isShotgun() && precursorMz>=0 ){
                  sc = new CgScan(0);
                  sc.Num = num;
                  sc.MsLevel = 1;
                  sc.RetentionTime = retentionTime;
                  sc.BasePeakMz = basePeakMz;
                  sc.BasePeakIntensity = basePeakIntensity;
                  sc.TotIonCurrent = totIonCurrent;
                  sc.setPolarity(polarity);
                  if (adders_ != null && adders_.length>0){                  
                    Vector<CgScan> scans = new Vector<CgScan>();
                    Vector<Range> ranges = new Vector<Range>();
                    for (AddScan adder : adders_){
                      if (adders_.length==1) scans.add(sc);
                      else scans.add(new CgScan(sc));
                      ranges.add(new Range(adder.getLowerThreshold(),adder.getUpperThreshold()));
                    }
                    
                    ////XmlReadPeaks(scans,ranges,maxRange,peaksCount,false,foundMzBorders);
                    Hashtable<Integer,Vector<Float>> mzValues = new Hashtable<Integer,Vector<Float>>();
                    Hashtable<Integer,Vector<Float>> intensities = new Hashtable<Integer,Vector<Float>>();
                    float mz = Float.parseFloat(getPrecursorMzString(precursorMz_));
                    //TODO: I do not know from where to get the intensity, since the precursorIntensity is always zero -> ask Kim
                    float intensity = totIonCurrent;
                    for (int k=0;k!=scans.size();k++){
                    	mzValues.put(k, new Vector<Float>());
                    	intensities.put(k, new Vector<Float>());
                      //TODO: I am not sure whether I should keep this check
                      if (ranges.get(k).getStart()<=mz && mz<ranges.get(k).getStop()){  
                        mzValues.get(k).add(mz);
                        intensities.get(k).add(intensity);
                      }
                    }
                    for (int k=0;k!=scans.size();k++){
                      CgScan aSc = scans.get(k);
                      aSc.PeaksCount = mzValues.get(k).size();
                      aSc.Scan = new float[mzValues.get(k).size()][2];
                      for (i=0; i!=mzValues.get(k).size();i++){
                        aSc.Scan[i][0] = mzValues.get(k).get(i);
                        aSc.Scan[i][1] = intensities.get(k).get(i);
                      }
                      if (aSc.PeaksCount>0 && !foundMzBorders){
                        aSc.LowMz =  mz*0.999f;
                        aSc.HighMz = mz*1.001f;
                      }
                    }
                    
                    for (int j=0; j!=adders_.length; j++){
                      AddScan adder = adders_[j];
                      sc = scans.get(j);
                      int currentLowMz = Math.round(sc.LowMz*multiplicationFactorForInt_);
                      int currentHighMz = Math.round(sc.HighMz*multiplicationFactorForInt_);
                      if (currentLowMz<this.lowestMz_) this.lowestMz_ = currentLowMz ;
                      if (currentHighMz>this.highestMz_) this.highestMz_ = currentHighMz;
                      // =================================================
                      // Now we can inform our caller that we have a valid
                      // new scan!
                      // =================================================
                      adder.AddScan(sc);
                    }
                  } else
                    throw new CgException("No adder for Header and Scans defined.");
                }*/
                if (baseScans==null){
                  baseScans = new Vector<CgScan>();
                  scanRanges = new Vector<Range>();
                  for (AddScan adder : this.adders_){
                    if (adder.getLastBaseScan()!=null){
                      baseScans.add(adder.getLastBaseScan());
                      scanRanges.add(new Range(adder.getLowerThreshold(),adder.getUpperThreshold()));
                    }
                  }
                }
                if (baseScans.size()>0){
                  if (parseMsMs_) {
                    MsMsScan msmsSc = new MsMsScan(peaksCount, num, msLevel, retentionTime,
                        lowMz, highMz, basePeakMz, basePeakIntensity,
                        totIonCurrent, getPrecursorMzString(precursorMz_), precursorIntensity,
                        polarity);
                    Vector<CgScan> qualifiedBaseScans = new Vector<CgScan>();
                    Vector<Range> qualifiedRanges = new Vector<Range>();
                    Vector<CgScan> ms2Scans = new Vector<CgScan>();
                    for (int j=0; j!=baseScans.size(); j++){
                      Range range = scanRanges.get(j);
                      float precMz = msmsSc.getMs1PrecursorMz();//Float.valueOf(precursorMz);
                      if (precMz<range.getStart() || range.getStop()<=precMz) continue;
                      //if (!range.insideRange(Float.valueOf(precursorMz))) continue;
                      qualifiedBaseScans.add(baseScans.get(j));
                      qualifiedRanges.add(range);
                      ms2Scans.add(new MsMsScan(msmsSc));
                    }
                    
                    if (qualifiedBaseScans.size()>0){
                      //TODO: in the msms subscans should not be any ranges
                      XmlReadPeaks(ms2Scans,qualifiedRanges,maxRange,peaksCount,true,foundMzBorders);
                      for (int j=0; j!=qualifiedBaseScans.size(); j++){
                        qualifiedBaseScans.get(j).AddSubscan(ms2Scans.get(j));
                      }
                      if (msLevel > myHeader_.highestMSLevel) myHeader_.highestMSLevel=msLevel;
                    }

                  } else {
                    for (CgScan scBase : baseScans) scBase.AddSubscanNumber(num);
                  }
                } else
                  throw new CgException("No base scan for subscan.");
              }
            } else if (parseMsMs_ && rdr_.getLocalName().equalsIgnoreCase("precursorMz")) {
              int attributeCount = rdr_.getAttributeCount();
              for (i = 0; i < attributeCount; i++) {
                if (rdr_.getAttributeLocalName(i) == "precursorIntensity") {
                  precursorIntensity = Float.parseFloat(rdr_.getAttributeValue(i));
                  try {
                    rdr_.next();
                    String childNode = rdr_.getText().trim();
                    if (childNode != null) {
                    //TODO: the following lines are for QQQ PIS/NLS - activate when necessary
/****                      precursorMz = */Float.parseFloat(childNode);
                      precursorMz_.add(childNode);
                      break;
                    }
                  }
                  catch (Exception ex) {
                    throw new CgException(ex.getMessage());
                  }
                }
              }
            } else if (rdr_.getLocalName().equalsIgnoreCase("scan")) {
              if (baseScans != null && baseScans.size()>0){
                XmlReadScan(baseScans,scanRanges,maxRange);
              }else{
                baseScans = new Vector<CgScan>();
                scanRanges = new Vector<Range>();
                for (AddScan adder : this.adders_){
                  if (adder.getLastBaseScan()!=null){
                    baseScans.add(adder.getLastBaseScan());
                    scanRanges.add(new Range(adder.getLowerThreshold(),adder.getUpperThreshold()));
                  }
                }
                XmlReadScan(baseScans,scanRanges,maxRange);
              }
            }
            break;
          case XMLStreamReader.END_ELEMENT:
            if (rdr_.getLocalName().equalsIgnoreCase("scan"))
              return;
            break;
        }
        eventType = rdr_.next();
      } while (eventType != XMLStreamReader.END_DOCUMENT);
    }
    catch (Exception ex) {
      ex.printStackTrace();
      throw new CgException(ex.getMessage());
    }
  }
  

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
  private void XmlReadPeaks(Vector<CgScan> scans, Vector<Range> ranges, Range maxRange, int peaksCount, boolean msms, boolean foundMzBorders) throws CgException
  {
    int i, j;
    String s;
    CgBase64 cgb = new CgBase64();

    // =========================================================
    // Read the peaks - Attributes:
    // =========================================================
    String compressionType = null;
    int precision = -1;
    
    if (scans != null && scans.size()>0) {
      String byteOrder = "";
      String pairOrder = "";

      for (i = 0; i < rdr_.getAttributeCount(); i++) {
        if (rdr_.getAttributeLocalName(i) == "precision") {
          precision = Integer.parseInt(rdr_.getAttributeValue(i));
        }
        if (rdr_.getAttributeLocalName(i) == "byteOrder") {
          byteOrder = rdr_.getAttributeValue(i);
        } else if (rdr_.getAttributeLocalName(i) == "pairOrder") {
          pairOrder = rdr_.getAttributeValue(i);
        }  else if (rdr_.getAttributeLocalName(i) == "compressionType") {
          compressionType = rdr_.getAttributeValue(i);
        }
      }
      for (CgScan scan : scans){
        scan.Precision = precision;
        scan.ByteOrder = byteOrder;
        scan.PairOrder = pairOrder;
      }
    }
    float lowestMzValue = Float.MAX_VALUE;
    float highestMzValue = 0f;
    try {
      rdr_.next();
      if (rdr_.getEventType()==XMLStreamReader.END_ELEMENT) return;
      s = rdr_.getText().trim(); // In s we have a Base64 coded value array!
      if (scans == null || scans.size()==0 || s == null || s.equalsIgnoreCase("</peaks>"))
        return;

      // =================================================
      // Process the data, in case we have to store it. In
      // C#, we have to store the byte array into a memory
      // stream and later on read it out float by float by
      // a binary reader.
      // =================================================
      byte[] decoded = cgb.decode(s);
      if (compressionType!=null && compressionType.equalsIgnoreCase("zlib")){
        Inflater decompressor = new Inflater();
        decompressor.setInput(decoded);
        ByteArrayOutputStream bos = null;
        try {
            bos = new ByteArrayOutputStream(decoded.length);

            // Decompress the data
            byte[] buf = new byte[1024];
            while (!decompressor.finished()) {
                int count = decompressor.inflate(buf);
                bos.write(buf, 0, count);
            }

        } finally {
            try {
                bos.close();
            } catch (Exception nope) { /* This exception doesn't matter */ }
        }
        decoded = bos.toByteArray();
      }
      
      ByteBuffer byteBuf = ByteBuffer.wrap(decoded);
      Hashtable<Integer,Vector<Float>> mzValues = new Hashtable<Integer,Vector<Float>>();
      Hashtable<Integer,Vector<Float>> intensities = new Hashtable<Integer,Vector<Float>>();
      for (int k=0;k!=scans.size();k++){
        mzValues.put(k, new Vector<Float>());
        intensities.put(k, new Vector<Float>());
      }
      if (precision==64){
        double doubleArray[] = new double[decoded.length/8];
        DoubleBuffer doubleBuf = byteBuf.asDoubleBuffer();
        doubleBuf.get(doubleArray);
        j = 0;
        for (i = 0; i < peaksCount; i++) {
          float mzValue = (float)doubleArray[j++];
          float intensity = (float)doubleArray[j++];
          if (mzValue<lowestMzValue) lowestMzValue = mzValue;
          if (mzValue>highestMzValue) highestMzValue = mzValue;
          if (!msms && scans.size()>1 && (mzValue<maxRange.getStart() || maxRange.getStop()<=mzValue))continue;
          for (int k=0;k!=scans.size();k++){
            if (msms || (ranges.get(k).getStart()<=mzValue && mzValue<ranges.get(k).getStop())){
              mzValues.get(k).add(mzValue);
              intensities.get(k).add(intensity);
            }
          }
        }        
      }else{
        float floatArray[] = new float[decoded.length / 4];
        FloatBuffer floatBuf = byteBuf.asFloatBuffer();
        floatBuf.get(floatArray);
        j = 0;
        for (i = 0; i < peaksCount; i++) {
          float mzValue = floatArray[j++];
          float intensity = floatArray[j++];
          if (mzValue<lowestMzValue) lowestMzValue = mzValue;
          if (mzValue>highestMzValue) highestMzValue = mzValue;
          if (!msms && scans.size()>1 && (mzValue<maxRange.getStart() || maxRange.getStop()<=mzValue))continue;
          for (int k=0;k!=scans.size();k++){
            if (msms || (ranges.get(k).getStart()<=mzValue && mzValue<ranges.get(k).getStop())){  
              mzValues.get(k).add(mzValue);
              intensities.get(k).add(intensity);
            }
          }
        }
      }
      for (int k=0;k!=scans.size();k++){
        CgScan sc = scans.get(k);
        sc.PeaksCount = mzValues.get(k).size();
        sc.Scan = new float[mzValues.get(k).size()][2];
        for (i=0; i!=mzValues.get(k).size();i++){
          sc.Scan[i][0] = mzValues.get(k).get(i);
          sc.Scan[i][1] = intensities.get(k).get(i);
        }
        if (sc.PeaksCount>0 && !foundMzBorders){
          sc.LowMz = lowestMzValue;
          sc.HighMz = highestMzValue;
        }
      }
    }
    catch (Exception ex) {
      ex.printStackTrace();
      throw new CgException(ex.getMessage());
    }
  }
  
  /**
   * Function converts an Time Interval to a double value[s]. Do something
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
  private String getPrecursorMzString (Vector<String> precursorMzs){
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
