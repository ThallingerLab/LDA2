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

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Vector;

import javax.xml.stream.XMLStreamException;
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
 * StAX based mzML Reader.
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
  
  //the required XML tag, of which the content will be read. 
  private static final String TAG_RUN = "run";
  
  //accessed XML sub-tags of TAG_RUN
  private static final String TAG_SPECTRUM_LIST = "spectrumList"; 
  private static final String TAG_SPECTRUM = "spectrum";
  private static final String TAG_CV_PARAM = "cvParam";
  private static final String TAG_SCAN_LIST = "scanList";
  private static final String TAG_PRECURSOR_LIST = "precursorList";
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
  private static final String ENTRY_LOWEST_OBSERVED_MZ = "lowest observed m/z"; //TODO: consider replacing with 'scan window lower limit'
  private static final String ENTRY_HIGHEST_OBSERVED_MZ = "highest observed m/z"; //TODO: consider replacing with 'scan window upper limit'
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
    String parentFileName = "";

    try {
      while (reader_.hasNext()) 
      {
      	if (isStartElement(TAG_RUN))
      	{
      		parentFileName = getRequiredAttribute(reader_, ATTRIBUTE_ID);
      		reader_.next();
      	}
      	else if (isStartElement(TAG_SPECTRUM_LIST))
      	{
      		int scanCount = Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_COUNT));
          informAdders(readOnlyRequiredInfoForMultiThreading, parentFileName, scanCount);
          readSpectrumList(readOnlyRequiredInfoForMultiThreading);
      	}
      	else if (isEndElement(this.getTagRun()))
      	{
      		return;
      	}
      	else
      	{
      		reader_.next();
      	}
      }
    }
    catch (Exception ex) {
      throw new CgException(ex.getMessage());
    }
  }
  
  /**
   * This method reads the full structure of the XML element constant TAG_SPECTRUM_LIST into memory.
   * @param readOnlyRequiredInfoForMultiThreading
   * 					Reads only the information required for getting an idea about the 
   *          m/z ranges for splitting in iterations and/or threads,
   *          the retention time range, 
   *          highest msLevel, 
   *          as well as whether polarity switching is used
   * @throws CgException
   * 					All internal exceptions are mapped to the CgException type.
   */
  private void readSpectrumList(boolean readOnlyRequiredInfoForMultiThreading) throws CgException
  {
  	AddScanHelper helper = new AddScanHelper();
  	Precursor precursor = new Precursor();
  	try {
      while (reader_.hasNext()) 
      {
      	if (isStartElement(TAG_SPECTRUM))
      	{
      		if (readOnlyRequiredInfoForMultiThreading) 
      		{
            readOnlyRequiredInfoForMultiThreading();
          } 
      		else 
          {
            readSpectrum(helper, precursor);
          }
      	}
      	else if (isEndElement(TAG_SPECTRUM_LIST))
      	{
      		return;
      	}
      	else
      	{
      		reader_.next();
      	}
      }
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
  	SpectrumParams params = new SpectrumParams();
  	boolean spectrumTagFound = false;
    //the following line is for QQQ PIS/NLS and PRM
//    float precursor = -1;
    
    try {
      while (reader_.hasNext()) {
      	if (isStartElement(TAG_SPECTRUM))
      	{
      		params.setSpectrumIndex(Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_INDEX)) +1); //mzXML files start at 1, mzML at 0
      		params.setPeaksCount(Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_DEFAULT_ARRAY_LENGTH)));
      		spectrumTagFound = true;
          reader_.next();
      	}
      	else if (spectrumTagFound && isStartElement(TAG_CV_PARAM))
      	{
      		parseSpectrumCVParamsOverview(params); //TODO: put the reader next thing in there, so that this part stops when there's a new section (e.g. scanlist)...
      		reader_.next();
      	}
      	else if (spectrumTagFound && isStartElement(TAG_BINARY_DATA_ARRAY))
      	{
      		boolean foundMzBorders = params.isFoundMzBorders();
        	evaluateScan(params);
        	
          if (params.getMsLevel() == 1) {
            if (!foundMzBorders && (params.getPeaksCount()>0)){
              float[] maxima = getMaximaFromBinaryDataArray(params.getPeaksCount());
              if (!Arrays.equals(maxima, new float[]{Float.MAX_VALUE, 0f})) {
              	params.setLowMz(maxima[0]);
              	params.setHighMz(maxima[1]);
              	params.setLowMzFound(true);
              	params.setHighMzFound(true);
              }
            }
            setCurrentGlobalMaxima(params.getLowMz(), params.getHighMz());
          }
          reader_.next();
      	}
      	else if (isEndElement(TAG_SPECTRUM))
      	{
      		if (params.getPolarity() == CgDefines.POLARITY_NO) {
            throw new CgException(String.format("The polarity at spectrum index %s could not be found.", 
                params.getSpectrumIndex()));
          } else {
            setPolaritySwitching(checkPolaritySwitching(params.getPolarity()));  
          }
      		return;
      	}
      	else
      	{
      		reader_.next();
      	}
      }
    }
    catch (Exception ex) {
      throw new CgException(ex.getMessage());
    }
  }
  
  /**
   * Parses an attribute to return the corresponding integer value.
   * @param entry							The entry (name) this attribute belongs to.
   * @param attribute					The attribute to parse the value of.
   * @return
   * @throws IOException			If the value of the attribute is missing.
   */
  private int parseIntegerValue(String entry, String attribute) throws IOException
  {
  	String value = reader_.getAttributeValue(null, attribute);
    if (!value.equals("")) 
    {
    	return Integer.parseInt(value);
    }
    else
    {
    	throw new IOException(String.format("A value for a '%s' entry is missing!", entry));
    }
  }
  
  /**
   * Parses an attribute to return the corresponding float value.
   * @param entry							The entry (name) this attribute belongs to (for error messages).
   * @param attribute					The attribute to parse the value of.
   * @return
   * @throws IOException			If the value of the attribute is missing.
   */
  private float parseFloatValue(String entry, String attribute) throws IOException
  {
  	String value = reader_.getAttributeValue(null, attribute);
    if (!value.equals("")) 
    {
    	return Float.parseFloat(value);
    }
    else
    {
    	throw new IOException(String.format("A value for a '%s' entry is missing!", entry));
    }
  }
  
  /**
   * This method reads out the m/z maxima of a binary data array, for getting an idea
   * about the m/z ranges for splitting in iterations and/or threads
   * 
   * @param peaksCount
   * @return A float[] containing the lowest and highest m/z value
   * @throws CgException
   */
  protected float[] getMaximaFromBinaryDataArray(int peaksCount) throws CgException
  {
    float[] maxima = new float[]{Float.MAX_VALUE, 0f};
    if (peaksCount<1) return maxima;
    String precision = null;
    String compression = null;
    
    try {
      while (reader_.hasNext()) {
      	if (isStartElement(TAG_CV_PARAM))
      	{
      		switch (getRequiredAttribute(reader_, ATTRIBUTE_NAME)) 
      		{
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
            case ENTRY_INTENSITY_ARRAY:
              return maxima;
            default:
            	break;
          }
      		reader_.next();
      	}
      	else if (isStartElement(TAG_BINARY))
      	{
          if (reader_.next() != XMLStreamReader.CHARACTERS) 
            return maxima;
          
          String binaryString = reader_.getText().trim();
          float floatArray[];
          if (precision==ENTRY_64_BIT_FLOAT)
          {
            floatArray = decode64(binaryString, compression);                  
          }
          else
          {
            floatArray = decode32(binaryString, compression); 
          }
          maxima[0] = floatArray[0];
          maxima[1] = floatArray[floatArray.length-1];
          reader_.next();
      	}
      	else if (isEndElement(TAG_BINARY_DATA_ARRAY))
      	{
      		return maxima;
      	}
      	else
      	{
      		reader_.next();
      	}
      }
    }
    catch (Exception ex) {
      throw new CgException(ex.getMessage());
    }   
    return maxima;
  }
  
  /**
   * This method decodes a float64 encoded binary String.
   * 
   * @param binary 				binary String to decode
   * @param compression 	compression of the binary String
   * 
   * @return a float array containing the decoded information
   */
  private float[] decode64(String binary, String compression) 
  {
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
   * @param binary 				binary String to decode
   * @param compression 	ompression of the binary String
   * 
   * @return a float array containing the decoded information
   */
  private float[] decode32(String binary, String compression) 
  {
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
   * Decodes a binary String with the appropriate method depending on the encoding precision.
   * @param precision			precision of the binary String
   * @param binaryString	binary String to decode	
   * @param compression		compression of the binary String
   * @return a float array containing the decoded information
   */
  private float[] decodeBinaryArray(String precision, String binaryString, String compression)
  {
  	float binaryArray[] = {};
  	if (precision==ENTRY_64_BIT_FLOAT)
    {
      binaryArray = decode64(binaryString, compression);
    } 
    else 
    {
      binaryArray = decode32(binaryString, compression); 
    }
  	return binaryArray;
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
    myHeader_.highestMSLevel = 1;
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
   * Only supports time given in units of minutes and seconds. 
   * Needs to be extended should the need arise!
   * 
   * @param time 			String to convert
   * @param unit 			Unit of the given time
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
  
  private boolean isStartElement(String localName)
  {
  	return isElement(XMLStreamReader.START_ELEMENT, localName);
  }
  
  private boolean isEndElement(String localName)
  {
  	return isElement(XMLStreamReader.END_ELEMENT, localName);
  }
  
  private boolean isElement(int eventType, String localName)
  {
  	if (reader_.getEventType() == eventType &&
  			reader_.getLocalName().equalsIgnoreCase(localName))
  	{
  		return true;
  	}
  	return false;
  }
  
  
  /**
   * This method reads the full structure of the XML element constant TAG_SPECTRUM into memory.
   * 
   * @param helper					Object containing the CgScan that represents the level 1 scan to which this scan belongs
   * 												as well as the individually allowed m/z ranges of the AddScan interfaces.
   * @param precursor				Object containing relevant information about precursors at lower MS levels.
   * @throws CgException		All internal exceptions are mapped to the CgException type.
   */
  protected void readSpectrum(AddScanHelper helper, Precursor precursor) throws CgException
  {
  	SpectrumParams params = new SpectrumParams();
    CgScan sc = null;
    Vector<Float> scanStartTimes = new Vector<Float>();
    boolean spectrumTagFound = false;
    
    try {
      while (reader_.hasNext()) {
      	if (isStartElement(TAG_SPECTRUM))
      	{
      		parseSpectrum(params);
      		spectrumTagFound = true;
      	}
      	else if (spectrumTagFound && isStartElement(TAG_SCAN_LIST)) //TODO: I do not have any examples with more than one scan here, thus not tested for more than one
      	{
      		parseScanList(scanStartTimes, params);
      	}
      	else if (spectrumTagFound && isStartElement(TAG_PRECURSOR_LIST))
      	{
      		parsePrecursorList(precursor, params.getMsLevel());
      	}
      	else if (spectrumTagFound && isStartElement(TAG_BINARY_DATA_ARRAY_LIST))
      	{
      		evaluateScan(params);
      		
      		if (params.getMsLevel() == 1) 
      		{
      			readMS1BinaryDataArray(sc, params, precursor, scanStartTimes);
      		}
      		else
      		{
      		  //the following lines (until end of bracket) are for QQQ PIS/NLS and PRM data
            //this generates an artificial MS1 scans for shotgun MSn-only data
            if (!myHeader_.hasMS1Scans && LipidomicsConstants.isShotgun()>LipidomicsConstants.SHOTGUN_FALSE && precursor.getPrecursorMzFloat()>=0 )
            {
              //TODO: do shotgun stuff
            }
            readMSnBinaryDataArray(params, precursor, scanStartTimes, helper);
            reader_.next();
      		}
      	}
      	
      	else if (isEndElement(TAG_SPECTRUM))
      	{
      		if (params.getPolarity() == CgDefines.POLARITY_NO) 
          {
            throw new CgException(String.format("The polarity at spectrum index %s could not be found.", 
            		params.getSpectrumIndex()));
          } else 
          {
            setPolaritySwitching(checkPolaritySwitching(params.getPolarity()));  
          }
      		helper.init(adders_);
      		return;
      	}
      	else
      	{
      		reader_.next();
      	}
      }
    }
    catch (Exception ex) {
      throw new CgException(ex.getMessage());
    }
  }
  
  
  /**
   * This method parses the structure of the XML element constant TAG_BINARY_DATA_ARRAY_LIST, belonging to a scan of an MS level higher than 1. 
   * 
   * @param params						Object containing the parameters of the corresponding spectrum.
   * @param precursor					Object containing relevant information about precursors at lower MS levels.
   * @param scanStartTimes		Vector containing the start times of the scan(s) belonging to this spectrum.
   * @param helper						Object containing the CgScan that represents the level 1 scan to which this scan belongs
   * 													as well as the individually allowed m/z ranges of the AddScan interfaces.
   * @throws CgException			All internal exceptions are mapped to the CgException type.
   */
  private void readMSnBinaryDataArray(SpectrumParams params, Precursor precursor, Vector<Float> scanStartTimes, AddScanHelper helper) throws CgException
  {
  	helper.init(adders_);
  	
  	if (helper.getScans().size()>0){
      if (this.getParseMsMs()) {
        MsMsScan msmsSc = new MsMsScan(params.getPeaksCount(), params.getSpectrumIndex(), params.getMsLevel(), 
        		scanStartTimes.get(0), //TODO: test data never has more than one scan, adapt in case more than one scan at different retention times happens...
        		params.getLowMz(), params.getHighMz(), params.getBasePeakMz(), params.getBasePeakIntensity(),
        		params.getTotalIonCurrent(), getPrecursorMzString(precursor.getPrecursorMz()), 0f, //setting precursor intensity to 0f, as it is never used
        		params.getPolarity());
        Vector<CgScan> qualifiedBaseScans = new Vector<CgScan>();
        Vector<Range> qualifiedRanges = new Vector<Range>();
        Vector<CgScan> ms2Scans = new Vector<CgScan>();
        for (int j=0; j!=helper.getScans().size(); j++){
          Range range = helper.getRanges().get(j);
          float precMz = msmsSc.getMs1PrecursorMz();
          if (!range.insideRange(Float.valueOf(precMz))) continue;
          qualifiedBaseScans.add(helper.getScans().get(j));
          qualifiedRanges.add(range);
          ms2Scans.add(new MsMsScan(msmsSc));
        }
        
        if (qualifiedBaseScans.size()>0){
          //TODO: in the msms subscans should not be any ranges
          readPeaks(ms2Scans,qualifiedRanges,true,params);
          for (int j=0; j!=qualifiedBaseScans.size(); j++){
            qualifiedBaseScans.get(j).AddSubscan(ms2Scans.get(j));
          }
          if (params.getMsLevel() > myHeader_.highestMSLevel) myHeader_.highestMSLevel=params.getMsLevel();
        }

      } else {
        for (CgScan scBase : helper.getScans()) scBase.AddSubscanNumber(params.getSpectrumIndex());
      }
    } 
    else
    {
    	throw new CgException("No base scan for subscan.");
    }
  }
  
  
  /**
   * This method parses the structure of the XML element constant TAG_BINARY_DATA_ARRAY_LIST, belonging to a scan of MS level = 1. 
   * 
   * @param sc								CgScan Object combining all information from this spectrum.					
   * @param params						Object containing the parameters of the corresponding spectrum.
   * @param precursor					Object containing relevant information about precursors at lower MS levels.
   * @param scanStartTimes		Vector containing the start times of the scan(s) belonging to this spectrum.
   * @throws CgException			All internal exceptions are mapped to the CgException type.
   */
  private void readMS1BinaryDataArray(CgScan sc, SpectrumParams params, Precursor precursor, Vector<Float> scanStartTimes) throws CgException
  {
		sc = new CgScan(0, params.getSpectrumIndex()+1, params.getMsLevel(), 
				scanStartTimes.get(0), //TODO: test data never has more than one scan, adapt in case more than one scan at different retention times happens...
				params.getLowMz(), params.getHighMz(), params.getBasePeakMz(), params.getBasePeakIntensity(),
    		params.getTotalIonCurrent(), params.getPolarity());
		
		if (adders_ != null && adders_.length>0){
			AddScanHelper helper = new AddScanHelper();
			helper.init(adders_, sc);
      readPeaks(helper.getScans(),helper.getRanges(),false,params);
      
      for (int j=0; j!=adders_.length; j++){
        AddScan adder = adders_[j];
        sc = helper.getScans().get(j);
        setCurrentGlobalMaxima(sc.LowMz, sc.HighMz);
       
        // =================================================
        // Now we can inform our caller that we have a valid
        // new scan!
        // =================================================
        adder.AddScan(sc);
      }
    } 
		else
    {
    	throw new CgException("No adder for Header and Scans defined.");
    } 
  }
  
  
  /**
   * This method parses the top level information of a an XML element constant TAG_SPECTRUM.
   * 
   * @param params					Object to hold the parsed information.
   * @throws CgException		All internal exceptions are mapped to the CgException type.
   */
  private void parseSpectrum(SpectrumParams params) throws CgException
  {
  	params.setSpectrumIndex(Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_INDEX)) +1); //mzXML files start at 1, mzML at 0
		params.setPeaksCount(Integer.parseInt(getRequiredAttribute(reader_, ATTRIBUTE_DEFAULT_ARRAY_LENGTH)));
  	
  	try {
      while (reader_.hasNext()) {
      	if (reader_.next() == XMLStreamReader.START_ELEMENT)
      	{
        	if (reader_.getLocalName().equalsIgnoreCase(TAG_CV_PARAM)) 
          {
        		parseSpectrumCVParams(params);
          }
          else
          {
          	return;
          }
      	}
      }
  	} 
  	catch (XMLStreamException | IOException ex) 
  	{
  		throw new CgException(ex.getMessage());
  	}
  }
  
  /**
   * This method parses all TAG_CV_PARAM elements of an XML element constant TAG_SPECTRUM.
   * 
   * @param params					Object to hold the parsed information.
   * @throws IOException		Thrown if something went wrong parsing the file.
   */
  private void parseSpectrumCVParams(SpectrumParams params) throws IOException
  {
  	switch (getRequiredAttribute(reader_, ATTRIBUTE_NAME)) 
    {
      case ENTRY_MS_LEVEL:
      	params.setMsLevel(parseIntegerValue(ENTRY_MS_LEVEL, ATTRIBUTE_VALUE));
        break;
      case ENTRY_POSITIVE_SCAN:
      	params.setPolarity(CgDefines.POLARITY_POSITIVE);
        break;
      case ENTRY_NEGATIVE_SCAN:
      	params.setPolarity(CgDefines.POLARITY_NEGATIVE);
        break;
      case ENTRY_BASE_PEAK_MZ:
      	params.setBasePeakMz(parseFloatValue(ENTRY_BASE_PEAK_MZ, ATTRIBUTE_VALUE));
        break;
      case ENTRY_BASE_PEAK_INTENSITY:
      	params.setBasePeakIntensity(parseFloatValue(ENTRY_BASE_PEAK_INTENSITY, ATTRIBUTE_VALUE));
        break;
      case ENTRY_TOTAL_ION_CURRENT:
      	params.setTotalIonCurrent(parseFloatValue(ENTRY_TOTAL_ION_CURRENT, ATTRIBUTE_VALUE));
        break;
      case ENTRY_LOWEST_OBSERVED_MZ:
      	params.setLowMz(parseFloatValue(ENTRY_LOWEST_OBSERVED_MZ, ATTRIBUTE_VALUE));
      	params.setLowMzFound(true);
        break;
      case ENTRY_HIGHEST_OBSERVED_MZ:
      	params.setHighMz(parseFloatValue(ENTRY_HIGHEST_OBSERVED_MZ, ATTRIBUTE_VALUE));
      	params.setHighMzFound(true);
        break;
      default:
      	break;
    }
  }
  
  /**
   * This method parses all TAG_CV_PARAM elements of an XML element constant TAG_SPECTRUM,
   * which are relevant for multithreading.
   * 
   * @param params					Object to hold the parsed information.
   * @throws IOException		Thrown if something went wrong parsing the file.
   */
  private void parseSpectrumCVParamsOverview(SpectrumParams params) throws IOException
  {
  	switch (getRequiredAttribute(reader_, ATTRIBUTE_NAME)) 
    {
      case ENTRY_MS_LEVEL:
      	params.setMsLevel(parseIntegerValue(ENTRY_MS_LEVEL, ATTRIBUTE_VALUE));
        break;
      case ENTRY_POSITIVE_SCAN:
      	params.setPolarity(CgDefines.POLARITY_POSITIVE);
        break;
      case ENTRY_NEGATIVE_SCAN:
      	params.setPolarity(CgDefines.POLARITY_NEGATIVE);
        break;
      case ENTRY_LOWEST_OBSERVED_MZ:
      	params.setLowMz(parseFloatValue(ENTRY_LOWEST_OBSERVED_MZ, ATTRIBUTE_VALUE));
      	params.setLowMzFound(true);
        break;
      case ENTRY_HIGHEST_OBSERVED_MZ:
      	params.setHighMz(parseFloatValue(ENTRY_HIGHEST_OBSERVED_MZ, ATTRIBUTE_VALUE));
      	params.setHighMzFound(true);
        break;
      case ENTRY_SCAN_START_TIME:
      	if (params.getSpectrumIndex() == 0) 
        {
          myHeader_.StartTime = convertTimeFormat(
              reader_.getAttributeValue(null, ATTRIBUTE_VALUE), reader_.getAttributeValue(null, ATTRIBUTE_UNIT_NAME));
        } 
        else if (params.getSpectrumIndex() == (myHeader_.ScanCount -1)) 
        {
          myHeader_.EndTime = convertTimeFormat(
              reader_.getAttributeValue(null, ATTRIBUTE_VALUE), reader_.getAttributeValue(null, ATTRIBUTE_UNIT_NAME));
        }
        break;
      default:
      	break;
    }
  }
  
  
  /**
   * This method parses the structure of the XML element constant TAG_SCAN_LIST. 
   * 
   * @param scanStartTimes			Vector containing the start times of the scan(s) belonging to this spectrum.
   * @param params							Object containing the parameters of the corresponding spectrum.
   * @throws CgException				All internal exceptions are mapped to the CgException type.
   */
  private void parseScanList(Vector<Float> scanStartTimes, SpectrumParams params) throws CgException
  {
  	try {
      while (reader_.hasNext()) 
      {
      	if (isStartElement(TAG_CV_PARAM) && getRequiredAttribute(reader_, ATTRIBUTE_NAME).equals(ENTRY_SCAN_START_TIME))
      	{
      		float scanStartTime = convertTimeFormat(
              reader_.getAttributeValue(null, ATTRIBUTE_VALUE), reader_.getAttributeValue(null, ATTRIBUTE_UNIT_NAME));
          scanStartTimes.add(scanStartTime);
          if (params.getSpectrumIndex() == 0) {
            myHeader_.StartTime = scanStartTime;
          } else if (params.getSpectrumIndex() == (myHeader_.ScanCount -1)) {
            myHeader_.EndTime = scanStartTime;
          }
          reader_.next();
      	}
      	else if (isEndElement(TAG_SCAN_LIST))
      	{
      		return;
      	}
      	else
      	{
      		reader_.next();
      	}
      }
  	} 
  	catch (XMLStreamException ex) 
  	{
  		throw new CgException(ex.getMessage());
  	}
  }
  
  
  /**
   * This method parses the structure of the XML element constant TAG_PRECURSOR_LIST. 
   * 
   * @param precursor				Object holding all precursor information.
   * @param msLevel					The MS level of the current spectrum.
   * @throws CgException		All internal exceptions are mapped to the CgException type.
   */
  private void parsePrecursorList(Precursor precursor, int msLevel) throws CgException
  {
  	if (msLevel<3)
  	{
  		precursor.setPrecursorMz(new Vector<String>());
    } 
  	else if (precursor.getLastMsLevel()>=msLevel)
  	{
      for (int i=(precursor.getLastMsLevel()+1); i>msLevel; i--){
        precursor.removeLastPrecursorMz();
      }
    }
    precursor.setLastMsLevel(msLevel);
  	
  	try {
  		while (reader_.hasNext()) {
  			if (isStartElement(TAG_CV_PARAM)) 
        {
      		switch (getRequiredAttribute(reader_, ATTRIBUTE_NAME)) 
          {
          	case ENTRY_SELECTED_ION_MZ:
          		precursor.addPrecursorMz(reader_.getAttributeValue(null, ATTRIBUTE_VALUE));
          		//the following line is for QQQ PIS/NLS and PRM
          		precursor.setPrecursorMzFloat(parseFloatValue(ENTRY_SELECTED_ION_MZ, ATTRIBUTE_VALUE));
          		reader_.next();
              break;
            default:
            	reader_.next();
            	break;
          }
        }
  			else if (isEndElement(TAG_PRECURSOR_LIST))
  			{
  				return;
  			}
  			else
  			{
  				reader_.next();
  			}
  		}
  	}
  	catch (XMLStreamException | IOException ex)
  	{
  		throw new CgException(ex.getMessage());
  	}
  }
  
  
  /**
   * This method evaluates the parameters of a given spectrum to set relevant global variables.
   * 
   * @param params			Object containing the parameters of the corresponding spectrum.
   */
  private void evaluateScan(SpectrumParams params)
  {
  	if (params.getMsLevel()==1)
  	{
  		if (!params.isFoundMzBorders() && !myHeader_.hasMS1Scans)
  		{
  			this.setLowestMz(ONE_MILLION*CgDefines.mzMultiplicationFactorForInt);
        this.setHighestMz(0); 
  		}
  		myHeader_.hasMS1Scans = true;
  	}

    if (params.getPolarity()!=CgDefines.POLARITY_NO)
    {
      if (currentPolarity_==CgDefines.POLARITY_NO)
      {
      	currentPolarity_ = params.getPolarity();
      }
      // if the polarity is suddenly different, polarity switching is used
      else if (currentPolarity_!=params.getPolarity())
      {
      	this.setPolaritySwitching(true);
      } 
    }
  }
  
  /**
   * This method reads the binary data arrays of a single scan and stores them in a CgScan object.
   * 
   * @param scans						Pass the CgScan object that represents the level 1 scan to which
   *          							this scan belongs to.
   * @param ranges					The m/z range restrictions that apply for every CgScan object
   * @param msms						True if the relevant scan is an MSn scan
   * @param params					Object containing the parameters of the corresponding spectrum.
   * @throws CgException		All internal exceptions are mapped to the CgException type.
   */
  protected void readPeaks(Vector<CgScan> scans, Vector<Range> ranges, boolean msms, SpectrumParams params) throws CgException
  {  
    float mzValueArray[] = {};
    float intensityValueArray[] = {};
    boolean mzValueArrayFound = false;
    boolean intensityValueArrayFound = false;
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
      while (reader_.hasNext()) {
      	if (isStartElement(TAG_CV_PARAM))
      	{
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
            default:
            	break;
          }
      	}
      	else if (isStartElement(TAG_BINARY) && valueType == ENTRY_MZ_ARRAY)
      	{
          if (reader_.next() != XMLStreamReader.CHARACTERS) return;
          
          mzValueArray = decodeBinaryArray(precision, reader_.getText().trim(), compression);
          mzValueArrayFound = true;
      	}
      	else if (isStartElement(TAG_BINARY) && valueType == ENTRY_INTENSITY_ARRAY)
      	{
          if (reader_.next() != XMLStreamReader.CHARACTERS) return;
          
          intensityValueArray = decodeBinaryArray(precision, reader_.getText().trim(), compression);
          intensityValueArrayFound = true;
      	}
      	else if (isEndElement(TAG_BINARY_DATA_ARRAY_LIST))
      	{
      		if (mzValueArrayFound&&intensityValueArrayFound) 
      		{
            for (int i = 0; i < params.getPeaksCount(); i++) 
            {
              float mzValue = mzValueArray[i];
              float intensityValue = intensityValueArray[i];
              if (!msms && scans.size()>1 && (mzValue<getMaxRange().getStart() || getMaxRange().getStop()<=mzValue)) continue;
              for (int k=0;k!=scans.size();k++)
              {
                if (msms || (ranges.get(k).getStart()<=mzValue && mzValue<ranges.get(k).getStop()))
                {  
                  mzValueHash.get(k).add(mzValue);
                  intensityValueHash.get(k).add(intensityValue);
                }
              }
            }
            for (int k=0;k!=scans.size();k++)
            {
              CgScan scan = scans.get(k);
              scan.PeaksCount = mzValueHash.get(k).size();
              scan.Scan = new float[scan.PeaksCount][2];
              for (int i=0; i!=scan.PeaksCount;i++)
              {
                scan.Scan[i][0] = mzValueHash.get(k).get(i);
                scan.Scan[i][1] = intensityValueHash.get(k).get(i);
              }
              if (scan.PeaksCount>0 && !params.isFoundMzBorders())
              {
                scan.LowMz = mzValueArray[0];
                scan.HighMz = mzValueArray[mzValueArray.length-1];
              }
            }
          }
          return;
      	}
        reader_.next();
      }
    }
    catch (Exception ex) {
      throw new CgException(ex.getMessage());
    } 
  }
  
  
  /**
   * Private helper class to encapsulate the parameters of a spectrum.
   * 
   * @author Leonida M. Lamp
   *
   */
  private class SpectrumParams
  {
		int spectrumIndex_ = -1;
    int peaksCount_ = 0;
  	int msLevel_ = 0;
  	int polarity_ = CgDefines.POLARITY_NO;
  	float basePeakMz_ = 0f;
  	float basePeakIntensity_ = 0f;
  	float totalIonCurrent_ = 0f;
  	float lowMz_ = Float.MAX_VALUE;
    float highMz_ = 0f;
    boolean lowMzFound_ = false;
    boolean highMzFound_ = false;
    
    public int getSpectrumIndex()
		{
			return spectrumIndex_;
		}

		public void setSpectrumIndex(int spectrumIndex)
		{
			this.spectrumIndex_ = spectrumIndex;
		}

		public int getPeaksCount()
		{
			return peaksCount_;
		}

		public void setPeaksCount(int peaksCount)
		{
			this.peaksCount_ = peaksCount;
		}

		public int getMsLevel()
		{
			return msLevel_;
		}

		public void setMsLevel(int msLevel)
		{
			this.msLevel_ = msLevel;
		}

		public int getPolarity()
		{
			return polarity_;
		}

		public void setPolarity(int polarity)
		{
			this.polarity_ = polarity;
		}

		public float getBasePeakMz()
		{
			return basePeakMz_;
		}

		public void setBasePeakMz(float basePeakMz)
		{
			this.basePeakMz_ = basePeakMz;
		}

		public float getBasePeakIntensity()
		{
			return basePeakIntensity_;
		}

		public void setBasePeakIntensity(float basePeakIntensity)
		{
			this.basePeakIntensity_ = basePeakIntensity;
		}

		public float getTotalIonCurrent()
		{
			return totalIonCurrent_;
		}

		public void setTotalIonCurrent(float totalIonCurrent)
		{
			this.totalIonCurrent_ = totalIonCurrent;
		}

		public float getLowMz()
		{
			return lowMz_;
		}

		public void setLowMz(float lowMz)
		{
			this.lowMz_ = lowMz;
		}

		public float getHighMz()
		{
			return highMz_;
		}

		public void setHighMz(float highMz)
		{
			this.highMz_ = highMz;
		}

		public void setLowMzFound(boolean lowMzFound)
		{
			this.lowMzFound_ = lowMzFound;
		}

		public void setHighMzFound(boolean highMzFound)
		{
			this.highMzFound_ = highMzFound;
		}
		
		public boolean isFoundMzBorders()
		{
			return this.highMzFound_ && this.lowMzFound_;
		}
  }
  
  /**
   * Private helper class to encapsulate CgScan objects and m/z range restrictions applying to them.
   * 
   * @author Leonida M. Lamp
   *
   */
  private class AddScanHelper
  {
  	Vector<CgScan> scans_ = null;
    Vector<Range> ranges_ = null;
    
    public void init(AddScan[] adders)
    {
    	if (scans_==null){
        scans_ = new Vector<CgScan>();
        ranges_ = new Vector<Range>();
        for (AddScan adder : adders){
          if (adder.getLastBaseScan()!=null){
            scans_.add(adder.getLastBaseScan());
            ranges_.add(new Range(adder.getLowerThreshold(),adder.getUpperThreshold()));
          }
        }
      }
    }
    
    public void init(AddScan[] adders, CgScan sc)
    {
    	scans_ = new Vector<CgScan>();
      ranges_ = new Vector<Range>();
      for (AddScan adder : adders){
        if (adders.length==1) scans_.add(sc);
        else scans_.add(new CgScan(sc));
        ranges_.add(new Range(adder.getLowerThreshold(),adder.getUpperThreshold()));
      }
    }

		public Vector<CgScan> getScans()
		{
			return scans_;
		}

		public Vector<Range> getRanges()
		{
			return ranges_;
		}
  }
  
  
  /**
   * Private helper class encapsulating all relevant precursor information.
   * 
   * @author Leonida M. Lamp
   *
   */
  private class Precursor
  {
  	Vector<String> precursorMz_ = new Vector<String>();
  	private int lastMsLevel_ = 0;
  	//the following line is for QQQ PIS/NLS and PRM data
  	float precursorMzFloat_ = -1f;
  	
		public Vector<String> getPrecursorMz()
		{
			return precursorMz_;
		}
		public void addPrecursorMz(String precursorMz)
		{
			this.precursorMz_.add(precursorMz);
		}
		public void setPrecursorMz(Vector<String> precursorMz)
		{
			this.precursorMz_ = precursorMz;
		}
		public void removeLastPrecursorMz()
		{
			precursorMz_.remove(precursorMz_.size()-1);
		}
		//the following method is for QQQ PIS/NLS and PRM data >>> not tested with mzML files
		public float getPrecursorMzFloat()
		{
			return precursorMzFloat_;
		}
		//the following method is for QQQ PIS/NLS and PRM data >>> not tested with mzML files
		public void setPrecursorMzFloat(float precursorMzFloat)
		{
			this.precursorMzFloat_ = precursorMzFloat;
		}
		public int getLastMsLevel()
		{
			return lastMsLevel_;
		}
		public void setLastMsLevel(int lastMsLevel)
		{
			this.lastMsLevel_ = lastMsLevel;
		}
		
  }
}