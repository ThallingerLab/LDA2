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

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.util.Hashtable;
import java.util.Vector;

import javax.xml.stream.XMLStreamReader;

import at.tugraz.genome.lda.LipidomicsConstants;
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
 * ATTENTION: Conversion from the raw data format has to be done with MsConvert, 
 * version number: 3.0.4976, as provided with LDA. Currently (August 2021), the official
 * release of MsConvert has a bug concerning the precursor masses for MS2 spectra.
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 */
public class MzXMLReader extends AbstractXMLSpectraReader
{ 
  //the call back interface for the header information
  private CgScanHeader myHeader_;
  
  /** the precursor m/z values of the current spectrum*/
  protected Vector<String> precursorMz_;
  
  /** the MS-level of the previous scan*/
  private int lastMsLevel_;
  
  //the required XML tag, of which the content will be read. 
  private static final String TAG_RUN = "msRun";
  
  //accessible XML sub-tags of TAG_RUN
  private static final String TAG_PARENT_FILE = "parentFile"; 
  private static final String TAG_SCAN = "scan";
  private static final String TAG_PEAKS = "peaks";
  private static final String TAG_END_PEAKS = "</peaks>";
  private static final String TAG_PRECURSOR_MZ = "precursorMz";
   
  //accessible XML attributes
  private static final String ATTRIBUTE_FILE_NAME = "fileName"; 
  private static final String ATTRIBUTE_SCAN_COUNT = "scanCount";
  private static final String ATTRIBUTE_START_TIME = "startTime";
  private static final String ATTRIBUTE_END_TIME = "endTime";
  private static final String ATTRIBUTE_NUM = "num";  
  private static final String ATTRIBUTE_MS_LEVEL = "msLevel";
  private static final String ATTRIBUTE_LOW_MZ = "lowMz";
  private static final String ATTRIBUTE_HIGH_MZ = "highMz";
  private static final String ATTRIBUTE_PEAKS_COUNT = "peaksCount";  
  private static final String ATTRIBUTE_POLARITY = "polarity";
  private static final String ATTRIBUTE_PRECURSOR_INTENSITY = "precursorIntensity";
  private static final String ATTRIBUTE_PRECISION = "precision";
  private static final String ATTRIBUTE_COMPRESSION_TYPE = "compressionType";
  private static final String ATTRIBUTE_RETENTION_TIME = "retentionTime";  
  private static final String ATTRIBUTE_BASE_PEAK_MZ = "basePeakMz";
  private static final String ATTRIBUTE_BASE_PEAK_INTENSITY = "basePeakIntensity";
  private static final String ATTRIBUTE_TOTAL_ION_CURRENT = "totIonCurrent";
  private static final String ATTRIBUTE_BYTE_ORDER = "byteOrder";
  private static final String ATTRIBUTE_PAIR_ORDER = "pairOrder";  
 
//  //accessible XML attribute-entries
  private static final String ENTRY_PLUS = "+";
  private static final String ENTRY_MINUS = "-";
  private static final String ENTRY_ZLIB = "zlib";
  
  /**
   * Constructs a MzMLReader object with given information.
   * 
   * @param callbacks
   *          Pass an array of objects implementing this interface. The methods will be
   *          called whenever a msScan header or individual scans are generated
   *          during the reading process.
   * @param parseMsMS Shall the MS/MS spectra be parsed too
   */
  public MzXMLReader(AddScan[] callbacks, boolean parseMsMs)
  {
    super(callbacks, parseMsMs);
  }
  
  /**
   * Constructs a MzMLReader object with given information.
   * 
   * @param callbacks
   *          Pass an array of objects implementing this interface. The methods will be
   *          called whenever a msScan header or individual scans are generated
   *          during the reading process.
   * @param parseMsMs Shall the MS/MS spectra be parsed too
   * @param multiplicationFactorForInt The multiplication factor to be used, to create integer values out of the m/z float values
   */
  public MzXMLReader(AddScan[] callbacks, boolean parseMsMs, int multiplicationFactorForInt)
  {
    super(callbacks, parseMsMs, multiplicationFactorForInt);
  }
  
  /**
   * Returns the XML element constant, the content of which will be read. 
   * Everything outside of it will be ignored by the reader.
   * 
   * @return the XML element constant
   */
  public String getTagRun()
  {
    return TAG_RUN;
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
  protected void readMsRun(boolean readOnlyRequiredInfoForMultiThreading) throws CgException
  {
    int i;
    int eventType;
    lastMsLevel_ = 0;
    
    myHeader_ = new CgScanHeader();
    writeMsRunAttributes(myHeader_);
    createAdderHeaders(myHeader_, readOnlyRequiredInfoForMultiThreading);

    try {
      do {
        eventType = reader_.getEventType();
        switch (eventType) {
          case XMLStreamReader.START_ELEMENT:
            if (reader_.getLocalName().equalsIgnoreCase(TAG_PARENT_FILE)){
              for (i = 0; i < reader_.getAttributeCount(); i++) {
                if (reader_.getAttributeLocalName(i).equalsIgnoreCase(ATTRIBUTE_FILE_NAME)) {
                  for (AddScan adder : adders_) adder.addParentFileName(StringUtils.getFileNameWOSuffix(reader_.getAttributeValue(i)));
                }
              }
            }else if (reader_.getLocalName().equalsIgnoreCase(TAG_SCAN)){
              if (readOnlyRequiredInfoForMultiThreading)
                readOnlyRequiredInfoForMultiThreading();
              else  
                readScan(null,null);
            }  
            break;
          case XMLStreamReader.END_ELEMENT:
            if (reader_.getLocalName().equalsIgnoreCase(getTagRun()))
              return;
            break;
        }
        eventType = reader_.next();
      } while (eventType != XMLStreamReader.END_DOCUMENT);
    }
    catch (Exception ex) {
      ex.printStackTrace();
      throw new CgException(ex.getMessage());
    }
  }
  
  /**
   * This method writes the MsRun attributes into a CgScanHeader object
   * 
   * @param myHeader_ The call back interface for the header information
   */
  private void writeMsRunAttributes(CgScanHeader myHeader_)
  {
    int count;
    myHeader_.highestMSLevel = 1;
    
    for (count = 0; count < reader_.getAttributeCount(); count++) {
      if (reader_.getAttributeLocalName(count).equalsIgnoreCase(ATTRIBUTE_SCAN_COUNT)) {
        myHeader_.ScanCount = Integer.parseInt(reader_.getAttributeValue(count));
      } 
      else if (reader_.getAttributeLocalName(count) == ATTRIBUTE_START_TIME) {
        myHeader_.StartTime = convertTimeFormat(reader_.getAttributeValue(count));
      } 
      else if (reader_.getAttributeLocalName(count) == ATTRIBUTE_END_TIME) {
        myHeader_.EndTime = convertTimeFormat(reader_.getAttributeValue(count));
      }
    } 
  }
  
  /**
   * Very basic function to convert the time string in mzXML files to a float value
   * 
   * @param s
   *          String to parse
   *          
   * @return float representing a number of seconds
   */
  protected float convertTimeFormat(String s)
  {
    if (s.startsWith("PT"))
      s = s.substring(2);
    if (s.endsWith("S"))
      s = s.substring(0, s.length() - 1);
    s = s.replace('e', 'E');
    return Float.parseFloat(s);
  }
  
  /**
   * This method creates headers for the AddScan objects
   * 
   * @param myHeader_ The call back interface for the header information
   * @param readOnlyRequiredInfoForMultiThreading 
   * @throws CgException
   */
  private void createAdderHeaders(CgScanHeader myHeader_, boolean readOnlyRequiredInfoForMultiThreading) throws CgException
  {
    if (adders_ != null && adders_.length>0)
      for (AddScan adder : adders_) adder.AddHeader(myHeader_);
    else{
      if (!readOnlyRequiredInfoForMultiThreading)
        throw new CgException("No adder for Header and Scans defined.");
    }     
  }
  
  /**
   * This method reads only the information required for getting an idea about the 
   * m/z ranges for splitting in iterations and/or threads,
   * the retention time range, 
   * highest msLevel, 
   * as well as whether polarity switching is used
   *          
   * @throws CgException All internal exceptions are mapped to the CgException type.
   */
  protected void readOnlyRequiredInfoForMultiThreading() throws CgException
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
    //the following line is for QQQ PIS/NLS and PRM
    float precursor = -1;

    for (i = 0; i < reader_.getAttributeCount(); i++) {
      if (reader_.getAttributeLocalName(i) == ATTRIBUTE_NUM) {
        num = Integer.parseInt(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_MS_LEVEL) {
        msLevel = Integer.parseInt(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_LOW_MZ) {
        lowMz = Float.parseFloat(reader_.getAttributeValue(i));
        lowMzFound = true;
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_HIGH_MZ) {
        highMz = Float.parseFloat(reader_.getAttributeValue(i));
        highMzFound=true;
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_PEAKS_COUNT) {
        peaksCount = Integer.parseInt(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_POLARITY) {
        polarityString = reader_.getAttributeValue(i);
        if (polarityString.equalsIgnoreCase(ENTRY_PLUS))
          polarity = CgDefines.POLARITY_POSITIVE;
        else if (polarityString.equalsIgnoreCase(ENTRY_MINUS))
          polarity = CgDefines.POLARITY_NEGATIVE;
        else
          throw new CgException(String.format("The scan contains an unknown polarity %s at scan number: %s", polarityString, num));
      }
    }
    boolean foundMzBorders = false;
    if (lowMzFound&&highMzFound) foundMzBorders = true;
    else if (msLevel==1 && !myHeader_.hasMS1Scans){
      this.setLowestMz(ONE_MILLION*CgDefines.mzMultiplicationFactorForInt);
      this.setHighestMz(0);      
    }
    if (msLevel==1)
      myHeader_.hasMS1Scans = true;

    if (polarity!=CgDefines.POLARITY_NO){
      if (currentPolarity_==CgDefines.POLARITY_NO)
        currentPolarity_ = polarity;
      // if the polarity is suddenly different, polarity switching is used
      else if (currentPolarity_!=polarity)
        this.setPolaritySwitching(true);
    }
    
    try {
      eventType = reader_.next();
      do {
        switch (eventType) {
          case XMLStreamReader.START_ELEMENT:
            if (reader_.getLocalName() == TAG_PEAKS) {
              if (msLevel == 1) {
                if (!foundMzBorders && (peaksCount>0)){
                  float[] maxima = getMaximaFromBinaryDataArray();
                  lowMz = maxima[0];
                  highMz = maxima[1];
                }
                setCurrentGlobalMaxima(lowMz, highMz);
                // =================================================
                // Now we can inform our caller that we have a valid
                // new scan!
                // =================================================
                
                //the following lines (until end of bracket) are for QQQ PIS/NLS and PRM data
              } else if (!myHeader_.hasMS1Scans && LipidomicsConstants.isShotgun()>LipidomicsConstants.SHOTGUN_FALSE && precursor>=0){
                setCurrentGlobalMaxima(precursor*0.999f, precursor*1.001f);
              }
              if (this.getParseMsMs() && msLevel > myHeader_.highestMSLevel) myHeader_.highestMSLevel=msLevel;
            //the following lines (until end of bracket) are for QQQ PIS/NLS and PRM data
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_PRECURSOR_MZ) && !myHeader_.hasMS1Scans && LipidomicsConstants.isShotgun()>LipidomicsConstants.SHOTGUN_FALSE) {
              int attributeCount = reader_.getAttributeCount();
              for (i = 0; i < attributeCount; i++) {
                if (reader_.getAttributeLocalName(i) == ATTRIBUTE_PRECURSOR_INTENSITY) {
                  try {
                    reader_.next();
                    String childNode = reader_.getText().trim();
                    if (childNode != null) {
                      precursor = Float.parseFloat(childNode);
                      break;
                    }
                  }
                  catch (Exception ex) {
                    throw new CgException(ex.getMessage());
                  }
                }
              }
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_SCAN)) {
              readOnlyRequiredInfoForMultiThreading();
            }
            break;
          case XMLStreamReader.END_ELEMENT:
            if (reader_.getLocalName().equalsIgnoreCase(TAG_SCAN))
              return;
            break;
        }
        eventType = reader_.next();
      } while (eventType != XMLStreamReader.END_DOCUMENT);
    }
    catch (Exception ex) {
      throw new CgException(ex.getMessage());
    }
    
  }
  
  /**
   * This method reads out the m/z maxima of a binary data array, for getting an idea
   * about the m/z ranges for splitting in iterations and/or threads
   * 
   * @return A float[] containing the lowest and highest m/z value
   * @throws CgException All internal exceptions are mapped to the CgException type.
   */
  protected float[] getMaximaFromBinaryDataArray() throws CgException
  {
    float[] maxima = new float[]{Float.MAX_VALUE,0f};
    int precision = -1;
    String compressionType = null;
    CgBase64 cgb = new CgBase64();
    for (int i = 0; i < reader_.getAttributeCount(); i++) {
      if (reader_.getAttributeLocalName(i).equalsIgnoreCase(ATTRIBUTE_PRECISION)) {
        precision = Integer.parseInt(reader_.getAttributeValue(i));
      }  else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_COMPRESSION_TYPE) {
        compressionType = reader_.getAttributeValue(i);
      }
    }
    try {
      reader_.next();
      String s = reader_.getText().trim(); // In s we have a Base64 coded value array!
      if (s == null || s.equalsIgnoreCase(TAG_END_PEAKS))
        return maxima;

      // =================================================
      // Process the data, in case we have to store it. In
      // C#, we have to store the byte array into a memory
      // stream and later on read it out float by float by
      // a binary reader.
      // =================================================
      byte[] decoded = cgb.decode(s);
      if (compressionType!=null && compressionType.equalsIgnoreCase(ENTRY_ZLIB)){
        decompressZLIB(decoded);
      }
      
      ByteBuffer byteBuf = ByteBuffer.wrap(decoded);
      float lowestMzValue = Float.MAX_VALUE;
      float highestMzValue = 0f;
      if (precision==64){
        double doubleArray[] = new double[decoded.length/8];
        DoubleBuffer doubleBuf = byteBuf.asDoubleBuffer();
        doubleBuf.get(doubleArray);
        lowestMzValue = (float)doubleArray[0];
        highestMzValue = (float)doubleArray[doubleArray.length-2];
      }else{
        float floatArray[] = new float[decoded.length / 4];
        FloatBuffer floatBuf = byteBuf.asFloatBuffer();
        floatBuf.get(floatArray);
        lowestMzValue = floatArray[0];
        highestMzValue = floatArray[floatArray.length-2];
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
   * This method reads a scan. It calls itself recursively in case scans are
   * element of our scan. However, if the scan level is > 1, the scan's content
   * is skipped for performance reasons.
   * 
   * @param scBase1
   *          Pass the CgScan object that represents the level 1 scan to which
   *          this scan belongs to.
   * @param ranges1 The individually allowed m/z ranges of the AddScan interfaces
   * 
   * @throws CgException All internal exceptions are mapped to the CgException type.
   */
  protected void readScan(Vector<CgScan> scBase1, Vector<Range> ranges1) throws CgException
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
    //the following line is for QQQ PIS/NLS and PRM data
    float precursorMz = -1f;

    // =========================================================
    // First of all we read the attributes:
    // =========================================================

    for (i = 0; i < reader_.getAttributeCount(); i++) {
      if (reader_.getAttributeLocalName(i) == ATTRIBUTE_NUM) {
        num = Integer.parseInt(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_MS_LEVEL) {
        msLevel = Integer.parseInt(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_PEAKS_COUNT) {
        peaksCount = Integer.parseInt(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_RETENTION_TIME) {
        retentionTime = convertTimeFormat(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_LOW_MZ) {
        lowMz = Float.parseFloat(reader_.getAttributeValue(i));
        lowMzFound = true;
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_HIGH_MZ) {
        highMz = Float.parseFloat(reader_.getAttributeValue(i));
        highMzFound = true;
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_BASE_PEAK_MZ) {
        basePeakMz = Float.parseFloat(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_BASE_PEAK_INTENSITY) {
        basePeakIntensity = Float.parseFloat(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_TOTAL_ION_CURRENT) {
        totIonCurrent = Float.parseFloat(reader_.getAttributeValue(i));
      } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_POLARITY) {
        polarityString = reader_.getAttributeValue(i);
        if (polarityString.equalsIgnoreCase(ENTRY_PLUS))
          polarity = CgDefines.POLARITY_POSITIVE;
        else if (polarityString.equalsIgnoreCase(ENTRY_MINUS))
          polarity = CgDefines.POLARITY_NEGATIVE;
        else
          throw new CgException(String.format("The scan contains an unknown polarity %s at scan number: %s", polarityString, num));
      }
    }
    boolean foundMzBorders = false;
    if (lowMzFound&&highMzFound) foundMzBorders = true;
    else if (msLevel==1 && !myHeader_.hasMS1Scans){
      this.setLowestMz(ONE_MILLION*CgDefines.mzMultiplicationFactorForInt);
      this.setHighestMz(0);      
    }
    if (msLevel==1)
      myHeader_.hasMS1Scans = true;
    if (polarity!=CgDefines.POLARITY_NO){
      if (currentPolarity_==CgDefines.POLARITY_NO)
        currentPolarity_ = polarity;
      // if the polarity is suddenly different, polarity switching is used
      else if (currentPolarity_!=polarity)
        this.setPolaritySwitching(true);
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
      eventType = reader_.next();
      do {
        switch (eventType) {
          case XMLStreamReader.START_ELEMENT:
            if (reader_.getLocalName().equalsIgnoreCase(TAG_PEAKS)) {
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

                  readPeaks(scans,ranges,peaksCount,false,foundMzBorders);
                  for (int j=0; j!=adders_.length; j++){
                    AddScan adder = adders_[j];
                    sc = scans.get(j);
                    setCurrentGlobalMaxima(sc.LowMz, sc.HighMz);
                    
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
                //the following lines (until end of bracket) are for QQQ PIS/NLS and PRM data
                //this generates an artificial MS1 scans for shotgun MSn-only data
                if (!myHeader_.hasMS1Scans && LipidomicsConstants.isShotgun()>LipidomicsConstants.SHOTGUN_FALSE && precursorMz>=0 ){
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
                      setCurrentGlobalMaxima(sc.LowMz, sc.HighMz);
                      
                      // =================================================
                      // Now we can inform our caller that we have a valid
                      // new scan!
                      // =================================================
                      adder.AddScan(sc);
                    }
                  } else
                    throw new CgException("No adder for Header and Scans defined.");
                }
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
                  if (this.getParseMsMs()) {
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
                      readPeaks(ms2Scans,qualifiedRanges,peaksCount,true,foundMzBorders);
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
            } else if (this.getParseMsMs() && reader_.getLocalName().equalsIgnoreCase(TAG_PRECURSOR_MZ)) {
              int attributeCount = reader_.getAttributeCount();
              for (i = 0; i < attributeCount; i++) {
                if (reader_.getAttributeLocalName(i) == ATTRIBUTE_PRECURSOR_INTENSITY) {
                  precursorIntensity = Float.parseFloat(reader_.getAttributeValue(i));
                  try {
                    reader_.next();
                    String childNode = reader_.getText().trim();
                    if (childNode != null) {
                    //the following line is for QQQ PIS/NLS and PRM data
                      precursorMz = Float.parseFloat(childNode);
                      precursorMz_.add(childNode);
                      break;
                    }
                  }
                  catch (Exception ex) {
                    throw new CgException(ex.getMessage());
                  }
                }
              }
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_SCAN)) {
              if (baseScans != null && baseScans.size()>0){
                readScan(baseScans,scanRanges);
              }else{
                baseScans = new Vector<CgScan>();
                scanRanges = new Vector<Range>();
                for (AddScan adder : this.adders_){
                  if (adder.getLastBaseScan()!=null){
                    baseScans.add(adder.getLastBaseScan());
                    scanRanges.add(new Range(adder.getLowerThreshold(),adder.getUpperThreshold()));
                  }
                }
                readScan(baseScans,scanRanges);
              }
            }
            break;
          case XMLStreamReader.END_ELEMENT:
            if (reader_.getLocalName().equalsIgnoreCase(TAG_SCAN))
              return;
            break;
        }
        eventType = reader_.next();
      } while (eventType != XMLStreamReader.END_DOCUMENT);
    }
    catch (Exception ex) {
      ex.printStackTrace();
      throw new CgException(ex.getMessage());
    }
  }
  

  /**
   * This method reads the binary data arrays of a single scan and stores them in a CgScan object.
   * 
   * @param scans
   *          Pass the CgScan object that represents the level 1 scan to which
   *          this scan belongs to.
   * @param ranges
   *          The m/z range restrictions that apply for every CgScan object
   * @param peaksCount The number of peaks in the binary data array
   * @param msms True if it is a MSn scan
   * @param foundMzBorders True if the lowest and highest m/z values of the file have been found
   * 
   * @throws CgException All internal exceptions are mapped to the CgException type.
   */
  protected void readPeaks(Vector<CgScan> scans, Vector<Range> ranges, int peaksCount, boolean msms, boolean foundMzBorders) throws CgException
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

      for (i = 0; i < reader_.getAttributeCount(); i++) {
        if (reader_.getAttributeLocalName(i) == ATTRIBUTE_PRECISION) {
          precision = Integer.parseInt(reader_.getAttributeValue(i));
        }
        if (reader_.getAttributeLocalName(i) == ATTRIBUTE_BYTE_ORDER) {
          byteOrder = reader_.getAttributeValue(i);
        } else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_PAIR_ORDER) {
          pairOrder = reader_.getAttributeValue(i);
        }  else if (reader_.getAttributeLocalName(i) == ATTRIBUTE_COMPRESSION_TYPE) {
          compressionType = reader_.getAttributeValue(i);
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
      reader_.next();
      if (reader_.getEventType()==XMLStreamReader.END_ELEMENT) return;
      s = reader_.getText().trim(); // In s we have a Base64 coded value array!
      if (scans == null || scans.size()==0 || s == null || s.equalsIgnoreCase(TAG_END_PEAKS))
        return;

      // =================================================
      // Process the data, in case we have to store it. In
      // C#, we have to store the byte array into a memory
      // stream and later on read it out float by float by
      // a binary reader.
      // =================================================    
      byte[] decoded = cgb.decode(s);
      if (compressionType!=null && compressionType.equalsIgnoreCase(ENTRY_ZLIB)){
        decompressZLIB(decoded);
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
          if (!msms && scans.size()>1 && (mzValue<getMaxRange().getStart() || getMaxRange().getStop()<=mzValue))continue;
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
          if (!msms && scans.size()>1 && (mzValue<getMaxRange().getStart() || getMaxRange().getStop()<=mzValue))continue;
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

}
