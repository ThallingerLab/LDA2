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
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.util.Arrays;
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


/**
 * StAX based mzML Reader for the Package.
 * 
 * This class reads an mzML file with hierarchical scans. It uses two methods
 * defined in the AddScan interface to return created header and spectrum
 * information to the caller. 
 * 
 * ATTENTION: Conversion from the raw data format has to be done with MsConvert, 
 * version number: 3.0.4976, as provided with LDA. Currently (August 2021), the official
 * release of MsConvert has a bug concerning the precursor masses for MS2 spectra.
 * 
 * @author Leonida M. Lamp
 */
public class MzMLReader extends AbstractXMLSpectraReader
{ 
  //the call back interface for the header information
  private CgScanHeader myHeader_;
  
  /** the precursor m/z values of the current spectrum*/
  protected Vector<String> precursorMz_;
  
  /** the MS-level of the previous scan*/
  private int lastMsLevel_;
  
  //the required XML tag, of which the content will be read. 
  private static final String TAG_RUN = "run";
  
  //accessible XML sub-tags of TAG_RUN
  private static final String TAG_SPECTRUM_LIST = "spectrumList"; 
  private static final String TAG_CV_PARAM = "cvParam";
  private static final String TAG_SPECTRUM = "spectrum";
  private static final String TAG_BINARY_DATA_ARRAY_LIST = "binaryDataArrayList";
  private static final String TAG_BINARY_DATA_ARRAY = "binaryDataArray";
  private static final String TAG_BINARY = "binary";
    
  //accessible XML attributes
  private static final String ATTRIBUTE_ID = "id"; 
  private static final String ATTRIBUTE_COUNT = "count";
  private static final String ATTRIBUTE_INDEX = "index";
  private static final String ATTRIBUTE_DEFAULT_ARRAY_LENGTH = "defaultArrayLength";
  private static final String ATTRIBUTE_NAME = "name";  
  private static final String ATTRIBUTE_VALUE = "value";
  private static final String ATTRIBUTE_UNIT_NAME = "unitName";
  
  //accessible XML attribute-entries
  private static final String ENTRY_MS_LEVEL = "ms level";
  private static final String ENTRY_LOWEST_OBSERVED_MZ = "lowest observed m/z";
  private static final String ENTRY_HIGHEST_OBSERVED_MZ = "highest observed m/z";
  private static final String ENTRY_POSITIVE_SCAN = "positive scan";
  private static final String ENTRY_NEGATIVE_SCAN = "negative scan";
  private static final String ENTRY_SCAN_START_TIME = "scan start time";
  private static final String ENTRY_64_BIT_FLOAT = "64-bit float";
  private static final String ENTRY_32_BIT_FLOAT = "32-bit float";
  private static final String ENTRY_NO_COMPRESSION = "no compression";
  private static final String ENTRY_ZLIB_COMPRESSION = "zlib compression";
  private static final String ENTRY_MZ_ARRAY = "m/z array";
  private static final String ENTRY_INTENSITY_ARRAY = "intensity array";
  private static final String ENTRY_BASE_PEAK_MZ = "base peak m/z";
  private static final String ENTRY_BASE_PEAK_INTENSITY = "base peak intensity";
  private static final String ENTRY_TOTAL_ION_CURRENT = "total ion current";
  private static final String ENTRY_MINUTE = "minute";
  private static final String ENTRY_SELECTED_ION_MZ = "selected ion m/z";
  private static final String ENTRY_PEAK_INTENSITY = "peak intensity";
  
  
   
  /**
   * Constructs a MzMLReader object with given information.
   * 
   * @param callbacks
   *          Pass an array of objects implementing the AddScan interface. The methods will be
   *          called whenever a msScan header or individual scans are generated during the reading
   *          process.
   * @param parseMsMS Shall the MS/MS spectra be parsed as well
   */
  public MzMLReader(AddScan[] callbacks, boolean parseMsMs)
  {
    super(callbacks, parseMsMs);
  }
  
  /**
   * Constructs a MzMLReader object with given information.
   * 
   * @param callbacks
   *          Pass an array of objects implementing the AddScan interface. The methods will be
   *          called whenever a msScan header or individual scans are generated during the reading 
   *          process.
   * @param parseMsMs Shall the MS/MS spectra be parsed as well
   * @param multiplicationFactorForInt The multiplication factor to be used in order to create integer values out of the m/z float values
   */
  public MzMLReader(AddScan[] callbacks, boolean parseMsMs, int multiplicationFactorForInt)
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
    int eventType;
    lastMsLevel_ = 0;
    String parentFileName = "";

    try {
      do {
        eventType = reader_.getEventType();
        switch (eventType) {
          
          case XMLStreamReader.START_ELEMENT:
            
            if (reader_.getLocalName().equalsIgnoreCase(TAG_RUN)){   
              parentFileName = getRequiredAttribute(reader_, ATTRIBUTE_ID);
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_SPECTRUM_LIST)) {
              int scanCount = Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_COUNT));
              informAdders(readOnlyRequiredInfoForMultiThreading, parentFileName, scanCount);
              if (readOnlyRequiredInfoForMultiThreading) {
                readOnlyRequiredInfoForMultiThreading();
              } else {
                readScan(null, null);
              }
            }
            break;
            
          case XMLStreamReader.END_ELEMENT:   
            
            if (reader_.getLocalName().equalsIgnoreCase(this.getTagRun()))
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
  
  //TODO: add stuff for shotgun (prm)
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
    int spectrumIndex = -1;
    int peaksCount = 0;
    int eventType;
    String value = null;
    int msLevel = 0;
    float lowMz = Float.MAX_VALUE;
    float highMz = 0f;
    boolean lowMzFound = false;
    boolean highMzFound = false;
    int polarity = CgDefines.POLARITY_NO;
    //the following line is for QQQ PIS/NLS and PRM
//    float precursor = -1;
    
    if (reader_.getLocalName().equalsIgnoreCase(TAG_SPECTRUM)) {
      spectrumIndex = Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_INDEX));
      peaksCount = Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_DEFAULT_ARRAY_LENGTH));
    }
    
    try {
      eventType = reader_.next();
      do {
        switch (eventType) {
          
          case XMLStreamReader.START_ELEMENT:
            
            if (reader_.getLocalName().equalsIgnoreCase(TAG_CV_PARAM)) {
              switch (getRequiredAttribute(reader_, ATTRIBUTE_NAME)) {
                case ENTRY_MS_LEVEL:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    msLevel = Integer.parseInt(value);
                    evaluateMsLevel(msLevel);
                  }
                  break;
                case ENTRY_LOWEST_OBSERVED_MZ:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    lowMz = Float.parseFloat(value);
                    lowMzFound = true;
                  }
                  break;
                case ENTRY_HIGHEST_OBSERVED_MZ:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    highMz = Float.parseFloat(value);
                    highMzFound = true;
                  }
                  break;
                case ENTRY_POSITIVE_SCAN:
                  polarity = CgDefines.POLARITY_POSITIVE;
                  break;
                case ENTRY_NEGATIVE_SCAN:
                  polarity = CgDefines.POLARITY_NEGATIVE; 
                  break;
                case ENTRY_SCAN_START_TIME:
                  if (spectrumIndex == 0) {
                    myHeader_.StartTime = convertTimeFormat(
                        reader_.getAttributeValue(null, ATTRIBUTE_VALUE), reader_.getAttributeValue(null, ATTRIBUTE_UNIT_NAME));
                  } else if (spectrumIndex == (myHeader_.ScanCount -1)) {
                    myHeader_.EndTime = convertTimeFormat(
                        reader_.getAttributeValue(null, ATTRIBUTE_VALUE), reader_.getAttributeValue(null, ATTRIBUTE_UNIT_NAME));
                  }
                  break;
              }
              
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_BINARY_DATA_ARRAY)) {
              if (msLevel == 1) {
                if (!(lowMzFound&&highMzFound) && (peaksCount>0)){
                  float[] maxima = getMaximaFromBinaryDataArray();
                  if (!Arrays.equals(maxima, new float[]{Float.MAX_VALUE, 0f})) {
                    lowMz = maxima[0];
                    highMz = maxima[1];
                    lowMzFound = highMzFound = true;
                  }
                }
                setCurrentGlobalMaxima(lowMz, highMz);
              }
              
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_SPECTRUM)) {
              readOnlyRequiredInfoForMultiThreading();
            }
            break;
            
          case XMLStreamReader.END_ELEMENT:
            
            if (reader_.getLocalName().equalsIgnoreCase(TAG_SPECTRUM)) {
              if (polarity == CgDefines.POLARITY_NO) {
                throw new CgException(String.format("The polarity at spectrum index %s could not be found.", 
                    spectrumIndex));
              } else {
                setPolaritySwitching(checkPolaritySwitching(polarity));  
              }
              return;   
              
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_SPECTRUM_LIST)){
              return;
            }
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
    float[] maxima = new float[]{Float.MAX_VALUE, 0f};
    float floatArray[];
    int eventType;
    String precision = null;
    String compression = null;
    
    try {
      eventType = reader_.next();
      
      do {
        switch (eventType) {
          
          case XMLStreamReader.START_ELEMENT:
            
            if (reader_.getLocalName().equalsIgnoreCase(TAG_CV_PARAM)) {
              switch (getRequiredAttribute(reader_, ATTRIBUTE_NAME)) {
                case ENTRY_64_BIT_FLOAT:
                  precision = ENTRY_64_BIT_FLOAT;
                  break;
                case ENTRY_32_BIT_FLOAT:
                  precision = ENTRY_32_BIT_FLOAT;
                  break;
                case ENTRY_NO_COMPRESSION:
                  compression = ENTRY_NO_COMPRESSION;
                  break;
                case ENTRY_ZLIB_COMPRESSION:
                  compression = ENTRY_ZLIB_COMPRESSION;
                  break;
                case ENTRY_MZ_ARRAY: 
                  break;
                case ENTRY_INTENSITY_ARRAY: // m/z maxima are stored in the m/z array
                  return maxima;
                  }
              
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_BINARY)) {
              eventType = reader_.next();
              if (eventType != XMLStreamReader.CHARACTERS) 
                return maxima;
              
              String binaryString = reader_.getText().trim();
 
              if (precision==ENTRY_64_BIT_FLOAT){
                floatArray = decode64(binaryString, compression);                  
              }else{
                floatArray = decode32(binaryString, compression); 
              }
              maxima[0] = floatArray[0];
              maxima[1] = floatArray[floatArray.length-1];
            }
            break;
          case XMLStreamReader.END_ELEMENT:
            if (reader_.getLocalName().equalsIgnoreCase(TAG_BINARY_DATA_ARRAY)) {
              return maxima;
            } 
            break;
        }
        eventType = reader_.next();
      } while (eventType != XMLStreamReader.END_DOCUMENT);
    }
    catch (Exception ex) {
      throw new CgException(ex.getMessage());
    }   
    return maxima;
  }
  
  /**
   * This method decodes a float64 encoded binary String.
   * 
   * @param binary The binary String to decode
   * @param compression The compression of the binary String
   * 
   * @return a float array containing the decoded information
   */
  private float[] decode64(String binary, String compression) {
    CgBase64 cgb = new CgBase64();    
    byte[] decoded = cgb.decode(binary);
    
    if (compression == ENTRY_ZLIB_COMPRESSION){
      decoded = decompressZLIB(decoded);
    }
    
    DoubleBuffer doubleBuf = ByteBuffer.wrap(decoded).order(ByteOrder.LITTLE_ENDIAN).asDoubleBuffer();
    double doubleArray[] = new double[decoded.length/8];
    doubleBuf.get(doubleArray);
    
    float[] floatArray = new float[doubleArray.length];
    for (int i = 0 ; i < doubleArray.length; i++)
        floatArray[i] = (float) doubleArray[i];
    
    return floatArray;    
  } 
  
  /**
   * This method decodes a float32 encoded binary String.
   * 
   * @param binary The binary String to decode
   * @param compression The compression of the binary String
   * 
   * @return a float array containing the decoded information
   */
  private float[] decode32(String binary, String compression) {
    CgBase64 cgb = new CgBase64();    
    byte[] decoded = cgb.decode(binary);
    
    if (compression == ENTRY_ZLIB_COMPRESSION){
      decoded = decompressZLIB(decoded);
    }
    
    FloatBuffer floatBuf = ByteBuffer.wrap(decoded).order(ByteOrder.LITTLE_ENDIAN).asFloatBuffer();
    float floatArray[] = new float[decoded.length/4];
    floatBuf.get(floatArray);
    
    return floatArray;    
  }  
  
  /**
   * This method checks if the polarity changes
   * 
   * @param polarity An integer representing the polarity
   */ 
  private boolean checkPolaritySwitching(int polarity) 
  {
    if (currentPolarity_==CgDefines.POLARITY_NO)
      currentPolarity_ = polarity;
    else if (currentPolarity_!=polarity)
      return true;
    return false;
  }
  
  /**
   * This method creates headers for the AddScan objects
   * 
   * @param readOnlyRequiredInfoForMultiThreading
   * @throws CgException
   */
  private void informAdders(boolean readOnlyRequiredInfoForMultiThreading, String parentFileName, int scanCount) throws CgException
  {
    myHeader_ = new CgScanHeader();
    myHeader_.ScanCount = scanCount;
    if (adders_ != null && adders_.length>0) 
      for (AddScan adder : adders_) {
        adder.AddHeader(myHeader_);
        adder.addParentFileName(parentFileName);
      }
    else{
      if (!readOnlyRequiredInfoForMultiThreading)
        throw new CgException("No adder for header and scans defined.");
    }     
  }
  
  /**
   * This method converts a String representing a time to a float.
   * Only supports time given in units of minutes and seconds so far, 
   * needs to be extended should the need arise!
   * 
   * @param time String to convert
   * @param unit Unit of the given time
   *          
   * @return float representing a number of seconds
   */
  protected float convertTimeFormat(String time, String unit)
  {
    float timeInSeconds;
    if (unit.equals(ENTRY_MINUTE))
      timeInSeconds = Float.parseFloat(time)*60;
    else
      timeInSeconds = Float.parseFloat(time);
    return timeInSeconds;
  }
  
  /**
   * This method sets the CgScanHeader information, lastMsLevel_ as well as 
   * the precursorMz_ Vector depending on the given msLevel.
   * 
   * @param msLevel MsLevel of the current scan.
   */
  private void evaluateMsLevel(int msLevel)
  {
    if (msLevel==1) {
      myHeader_.hasMS1Scans = true;
    }
    if (msLevel==1 || this.getParseMsMs()) {
      if (msLevel > myHeader_.highestMSLevel) 
        myHeader_.highestMSLevel=msLevel;
    }
    if (msLevel<3){
      precursorMz_ = new Vector<String>();
    } else if (lastMsLevel_>=msLevel){
      for (int i=(lastMsLevel_+1); i!=msLevel; i--){
        precursorMz_.remove(precursorMz_.size()-1);
      }
    }
    lastMsLevel_ = msLevel;
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
    Vector<CgScan> baseScans = scBase1;
    Vector<Range> scanRanges = ranges1;
    CgScan sc = null;
    String value = null;
    int spectrumIndex = -1;
    int peaksCount = 0;
    int eventType;
    int msLevel = 0;
    float lowMz = Float.MAX_VALUE;
    float highMz = 0f;
    boolean lowMzFound = false;
    boolean highMzFound = false;
    int polarity = CgDefines.POLARITY_NO;
    float scanStartTime = 0f;
    float basePeakMz = 0f;
    float basePeakIntensity = 0f;
    float totalIonCurrent = 0f;
    float precursorIntensity = 0f;
  //the following line is for QQQ PIS/NLS and PRM data
    float precursorMz = -1f;
    
    if (reader_.getLocalName().equalsIgnoreCase(TAG_SPECTRUM)) {
      spectrumIndex = Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_INDEX));
      peaksCount = Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_DEFAULT_ARRAY_LENGTH));
    }
    
    try {
      eventType = reader_.next();
      do {
        switch (eventType) {
          
          case XMLStreamReader.START_ELEMENT:
            
            if (reader_.getLocalName().equalsIgnoreCase(TAG_CV_PARAM)) {
              switch (getRequiredAttribute(reader_, ATTRIBUTE_NAME)) {
                case ENTRY_MS_LEVEL:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    msLevel = Integer.parseInt(value);
                    evaluateMsLevel(msLevel);
                  }
                  break;
                case ENTRY_LOWEST_OBSERVED_MZ:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    lowMz = Float.parseFloat(value);
                    lowMzFound = true;
                  }
                  break;
                case ENTRY_HIGHEST_OBSERVED_MZ:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    highMz = Float.parseFloat(value);
                    highMzFound = true;
                  }
                  break;
                case ENTRY_POSITIVE_SCAN:
                  polarity = CgDefines.POLARITY_POSITIVE;
                  break;
                case ENTRY_NEGATIVE_SCAN:
                  polarity = CgDefines.POLARITY_NEGATIVE; 
                  break;
                case ENTRY_SCAN_START_TIME:
                  scanStartTime = convertTimeFormat(
                      reader_.getAttributeValue(null, ATTRIBUTE_VALUE), reader_.getAttributeValue(null, ATTRIBUTE_UNIT_NAME));
                  if (spectrumIndex == 0) {
                    myHeader_.StartTime = scanStartTime;
                  } else if (spectrumIndex == (myHeader_.ScanCount -1)) {
                    myHeader_.EndTime = scanStartTime;
                  }
                  break;
                case ENTRY_BASE_PEAK_MZ:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    basePeakMz = Float.parseFloat(value);
                  }
                  break;
                case ENTRY_BASE_PEAK_INTENSITY:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    basePeakIntensity = Float.parseFloat(value);
                  }
                  break;
                case ENTRY_TOTAL_ION_CURRENT:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    totalIonCurrent = Float.parseFloat(value);
                  }
                  break;
                //TODO: mzML supports multiple precursor ions; handle that (precursorMz list?)!
                case ENTRY_SELECTED_ION_MZ:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    precursorMz = Float.parseFloat(value);
                    precursorMz_.add(value);
                  }
                  break; 
                case ENTRY_PEAK_INTENSITY:
                  value = reader_.getAttributeValue(null, ATTRIBUTE_VALUE);
                  if (!value.equals("")) {
                    precursorIntensity = Float.parseFloat(value);
                  }
                  break; 
              }
              
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_BINARY_DATA_ARRAY_LIST)) {
              if (msLevel == 1) {
                sc = new CgScan(0, spectrumIndex+1, msLevel, scanStartTime, lowMz, highMz, basePeakMz, 
                    basePeakIntensity, totalIonCurrent, polarity);
                
                if (adders_ != null && adders_.length>0){
                  Vector<CgScan> scans = new Vector<CgScan>();
                  Vector<Range> ranges = new Vector<Range>();
                  for (AddScan adder : adders_){
                    if (adders_.length==1) scans.add(sc);
                    else scans.add(new CgScan(sc));
                    ranges.add(new Range(adder.getLowerThreshold(),adder.getUpperThreshold()));
                  }
                  
                  if (scans != null && scans.size()!=0) {
                    readPeaks(scans,ranges,peaksCount,false,(lowMzFound&&highMzFound));
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
                  throw new CgException(
                      "No adder for Header and Scans defined.");
                
              } else {
                //the following lines (until end of bracket) are for QQQ PIS/NLS and PRM data
                //this generates an artificial MS1 scans for shotgun MSn-only data
                if (!myHeader_.hasMS1Scans && LipidomicsConstants.isShotgun()>LipidomicsConstants.SHOTGUN_FALSE && precursorMz>=0 ){
                  //TODO: do shotgun stuff
                }
                
                if (baseScans.size()>0){
                  if (this.getParseMsMs()) {
                    
                    MsMsScan msMsSc = new MsMsScan(peaksCount, spectrumIndex+1, msLevel, scanStartTime,
                        lowMz, highMz, basePeakMz, basePeakIntensity, totalIonCurrent, 
                        getPrecursorMzString(precursorMz_), precursorIntensity, polarity);
                    Vector<CgScan> qualifiedBaseScans = new Vector<CgScan>();
                    Vector<Range> qualifiedRanges = new Vector<Range>();
                    Vector<CgScan> ms2Scans = new Vector<CgScan>();
                    for (int j=0; j!=baseScans.size(); j++){
                      Range range = scanRanges.get(j);
                      float precMz = msMsSc.getMs1PrecursorMz();//Float.valueOf(precursorMz);
                      if (precMz<range.getStart() || range.getStop()<=precMz) continue;
                      //if (!range.insideRange(Float.valueOf(precursorMz))) continue;
                      qualifiedBaseScans.add(baseScans.get(j));
                      qualifiedRanges.add(range);
                      ms2Scans.add(new MsMsScan(msMsSc));
                    }
                    
                    if (qualifiedBaseScans.size()>0 && ms2Scans != null && ms2Scans.size()!=0){
                      //TODO: in the msms subscans should not be any ranges
                      readPeaks(ms2Scans,qualifiedRanges,peaksCount,true,(lowMzFound&&highMzFound));
                      
                      for (int j=0; j!=qualifiedBaseScans.size(); j++){
                        qualifiedBaseScans.get(j).AddSubscan(ms2Scans.get(j));
                      }
                    }
                  } else {
                    for (CgScan scBase : baseScans) scBase.AddSubscanNumber(spectrumIndex+1);
                  }
                } else
                  throw new CgException("No base scan for subscan.");     
              }
              
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_SPECTRUM)) {
              if (baseScans != null && baseScans.size()>0){
                readScan(baseScans,scanRanges);
              } else {
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
            
            if (reader_.getLocalName().equalsIgnoreCase(TAG_SPECTRUM)) {
              if (polarity == CgDefines.POLARITY_NO) {
                throw new CgException(String.format("The polarity at spectrum index %s could not be found.", 
                    spectrumIndex));
              } else {
                setPolaritySwitching(checkPolaritySwitching(polarity));  
              }
              return;  
              
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_SPECTRUM_LIST)){
              return;
            }
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
    float mzValueArray[] = {};
    float intensityValueArray[] = {};
    boolean mzValueArrayFound = false;
    boolean intensityValueArrayFound = false;
    int eventType;
    String precision = null;
    String compression = null;
    String valueType = "";
    Hashtable<Integer,Vector<Float>> mzValueHash = new Hashtable<Integer,Vector<Float>>();  
    Hashtable<Integer,Vector<Float>> intensityValueHash = new Hashtable<Integer,Vector<Float>>();
    
    for (int k=0;k!=scans.size();k++){
      mzValueHash.put(k, new Vector<Float>());
      intensityValueHash.put(k, new Vector<Float>());
    }
    
    try {
      eventType = reader_.next();
      
      do {
        switch (eventType) {
          
          case XMLStreamReader.START_ELEMENT:
            
            if (reader_.getLocalName().equalsIgnoreCase(TAG_CV_PARAM)) {
              switch (getRequiredAttribute(reader_, ATTRIBUTE_NAME)) {
                case ENTRY_64_BIT_FLOAT:
                  precision = ENTRY_64_BIT_FLOAT;
                  break;
                case ENTRY_32_BIT_FLOAT:
                  precision = ENTRY_32_BIT_FLOAT;
                  break;
                case ENTRY_NO_COMPRESSION:
                  compression = ENTRY_NO_COMPRESSION;
                  break;
                case ENTRY_ZLIB_COMPRESSION:
                  compression = ENTRY_ZLIB_COMPRESSION;
                  break;
                case ENTRY_MZ_ARRAY: 
                  valueType = ENTRY_MZ_ARRAY;
                  break;
                case ENTRY_INTENSITY_ARRAY:
                  valueType = ENTRY_INTENSITY_ARRAY;
                  break;
                  }
              
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_BINARY) && valueType == ENTRY_MZ_ARRAY) {
              eventType = reader_.next();
              if (eventType != XMLStreamReader.CHARACTERS) 
                return;
              
              String binaryString = reader_.getText().trim();
              if (precision==ENTRY_64_BIT_FLOAT){
                mzValueArray = decode64(binaryString, compression);
              } else {
                mzValueArray = decode32(binaryString, compression); 
              }
              mzValueArrayFound = true;
              
            } else if (reader_.getLocalName().equalsIgnoreCase(TAG_BINARY) && valueType == ENTRY_INTENSITY_ARRAY) {
              eventType = reader_.next();
              if (eventType != XMLStreamReader.CHARACTERS) 
                return;
              
              String binaryString = reader_.getText().trim();
              if (precision==ENTRY_64_BIT_FLOAT){
                intensityValueArray = decode64(binaryString, compression);
              } else {
                intensityValueArray = decode32(binaryString, compression); 
              }
              intensityValueArrayFound = true;
              
            }
            break;
            
          case XMLStreamReader.END_ELEMENT:
            
            if (reader_.getLocalName().equalsIgnoreCase(TAG_BINARY_DATA_ARRAY_LIST)) {
              if (mzValueArrayFound&&intensityValueArrayFound) {
                
                for (int i = 0; i < peaksCount; i++) {
                  float mzValue = mzValueArray[i];
                  float intensityValue = intensityValueArray[i];
                  if (!msms && scans.size()>1 && (mzValue<getMaxRange().getStart() || getMaxRange().getStop()<=mzValue))continue;
                  for (int k=0;k!=scans.size();k++){
                    if (msms || (ranges.get(k).getStart()<=mzValue && mzValue<ranges.get(k).getStop())){  
                      mzValueHash.get(k).add(mzValue);
                      intensityValueHash.get(k).add(intensityValue);
                    }
                  }
                }

                for (int k=0;k!=scans.size();k++){
                  CgScan scan = scans.get(k);
                  scan.PeaksCount = mzValueHash.get(k).size();
                  scan.Scan = new float[scan.PeaksCount][2];
                  for (int i=0; i!=scan.PeaksCount;i++){
                    scan.Scan[i][0] = mzValueHash.get(k).get(i);
                    scan.Scan[i][1] = intensityValueHash.get(k).get(i);
                  }
                  if (scan.PeaksCount>0 && !foundMzBorders){
                    scan.LowMz = mzValueArray[0];
                    scan.HighMz = mzValueArray[mzValueArray.length-1];
                  }
                }
                
              }
              return;
            } 
            break;
        }
        eventType = reader_.next();
      } while (eventType != XMLStreamReader.END_DOCUMENT);
    }
    catch (Exception ex) {
      throw new CgException(ex.getMessage());
    } 
  } 
  
  
}