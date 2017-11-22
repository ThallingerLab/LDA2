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

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.math.BigDecimal;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.dbutilities.Base64;
import at.tugraz.genome.maspectras.quantification.CgDefines;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgScan;
import at.tugraz.genome.maspectras.quantification.CgScanHeader;
import at.tugraz.genome.maspectras.quantification.MsMsScan;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.maspectras.utils.StringUtils;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class RawToChromThread extends Thread  implements AddScan
{

  /** indicating if the thread is finished*/
  private boolean finished_;
  /** the error string if something goes wrong*/
  private String errorString_;
  /** the directory where the chrom files have to be stored*/
  private String dir_;
  /** the lower m/z threshold for translating the data in integer format*/
  private Integer lowerThreshInt_;
  /** the upper m/z threshold for translating the data in integer format*/
  private Integer upperThreshInt_;
  /** the lower m/z threshold for data parsing*/
  private Float lowerThresholdTol_;
  /** the upper m/z threshold for data parsing*/
  private Float upperThresholdTol_;

  /** actual Scan Count (used as index) */
  private Hashtable<String,Integer> m_scanCount_;
  /** the header information*/
  private Hashtable<String,CgScanHeader> headerHash_;
  /** the individual scans*/
  private Hashtable<String,CgScan[]> scanHash_;
  /** the arraysize for reading m/z ranges in batch mode*/
  private int elementsForBatchCalculation_;
  
  /** allow reading of MSn information */
  private boolean msms_;
  /** is MSn information distributed in several files, such as in Waters MSE data*/
  private boolean msmsInSeveralFiles_;
  /** the multiplication factor to be used, to create integer values out of the m/z float values*/
  private int multiplicationFactorForInt_;
  /** the lowest resolution distance value by a given multiplicationFactorForInt_*/ 
  private int lowestResolution_;
  /** the lowest possible mz value*/
  private int lowestMz_;
  /** the highest possible mz value*/
  private int highestMz_;

  /** the highest MSn level in the data file*/
  private int highestMsLevel_;
  /** in my opinion, the msMsExclusionList was necessary because the scans were read overlapping - since the threaded version, this hashtable should be obsolete*/
  private Hashtable<String,Hashtable<Integer,Integer>> msMsExclusionList_;
  /** contains the amount of MSn scans*/
  private Hashtable<String,Hashtable<Integer,Integer>> numberOfMs2Scans_;
  
  /** the data file name (in case of merged data) that is currently read by the reader*/ 
  private String currentFileName_;
  
  /** if there are no merged MS runs, this value is used as key for the hashtables*/ 
  private final static String NO_NAME_PLACEHOLDER = "";
  
  /** the tolerance m/z value for overlapped reading*/
  private final static float ALLOWED_OVERLAP = 1f;

  
  /** the name of the chromatography file*/
  private Hashtable<String,String> chromFileName_;
  /** the name of the index file*/
  private Hashtable<String,String> indexFileName_;
  /** the name of the retention-time file*/
  private Hashtable<String,String> retentionTimeFileName_;

  /**
   * this constructor sets default values and inits required hash tables
   */
  private RawToChromThread(){
    finished_ = false;
    errorString_ = null;
    msMsExclusionList_ = new Hashtable<String,Hashtable<Integer,Integer>>();
    this.chromFileName_ = new Hashtable<String,String>();
    this.indexFileName_ = new Hashtable<String,String>();
    this.retentionTimeFileName_ = new Hashtable<String,String>();
    this.msms_ = false;
  }
  
  /**
   * 
   * @param msms allow reading of MSn information
   * @param multiplicationFactorForInt the multiplication factor to be used, to create integer values out of the m/z float values
   * @param lowestResolution the lowest resolution distance value by a given multiplicationFactorForInt
   */
  public RawToChromThread(boolean msms, int multiplicationFactorForInt, int lowestResolution, int highestMsLevel)
  {
    this();
    this.msms_ = msms;
    this.multiplicationFactorForInt_ = multiplicationFactorForInt;
    this.lowestResolution_ = lowestResolution;
    this.highestMsLevel_ = highestMsLevel;
  }
  
  /**
   * constructor creating a clone of the current object
   * @param otherThread the other RawToChromThread object to be cloned
   */
  public RawToChromThread(RawToChromThread otherThread){
    this();
    this.chromFileName_ = otherThread.chromFileName_;
    this.indexFileName_ = otherThread.indexFileName_;
    this.retentionTimeFileName_ = otherThread.retentionTimeFileName_;
    this.elementsForBatchCalculation_ = otherThread.elementsForBatchCalculation_;
    this.msms_ = otherThread.msms_;
    this.msmsInSeveralFiles_ = otherThread.msmsInSeveralFiles_;
    this.multiplicationFactorForInt_ = otherThread.multiplicationFactorForInt_;
    this.lowestResolution_ = otherThread.lowestResolution_;
    this.lowestMz_ = otherThread.lowestMz_;
    this.highestMz_ = otherThread.highestMz_;
    this.highestMsLevel_ = otherThread.highestMsLevel_;
  }
  
  public void run(){
    try{
      writeToChrom();
    } catch (Exception ex){
      ex.printStackTrace();
      errorString_ = ex.toString();
      
    }
    finished_ = true;
  }
  
  /**
   * sets the storage directory and the lower/upper m/z threshold for translating the data in integer format
   * @param dir the directory where the chrom files have to be stored
   * @param lowerThreshold the lower m/z threshold for translating the data in integer format
   * @param upperThreshold the upper m/z threshold for translating the data in integer format
   */
  public void setRequiredInformation(String dir, Integer lowerThreshold, Integer upperThreshold){
    this.dir_ = dir;
    this.lowerThreshInt_ = lowerThreshold;
    this.lowerThresholdTol_ = lowerThreshold.floatValue()/(float)this.multiplicationFactorForInt_-ALLOWED_OVERLAP;
    int lowerThreshIntTol = getCorrespondingIntValue(this.lowerThresholdTol_, multiplicationFactorForInt_);
    if (lowerThreshIntTol<this.lowestMz_){
      this.lowerThresholdTol_ = (float)this.lowestMz_/(float)this.multiplicationFactorForInt_;
    }
    this.upperThreshInt_ = upperThreshold;
    this.upperThresholdTol_ = upperThreshold.floatValue()/(float)this.multiplicationFactorForInt_+ALLOWED_OVERLAP;
    int upperThreshIntTol = getCorrespondingIntValue(this.upperThresholdTol_, multiplicationFactorForInt_);
    if (upperThreshIntTol>this.highestMz_) this.upperThresholdTol_ = (float)(this.highestMz_+1)/(float)this.multiplicationFactorForInt_; 
    this.initStoreHashes();
  }
  
  
  public void AddScan(CgScan sx) throws CgException
  {
    if (!scanHash_.containsKey(currentFileName_)) throw new CgException("m_scans Array not allocated");
    if (scanHash_.get(currentFileName_).length==m_scanCount_.get(currentFileName_)) return;
    scanHash_.get(currentFileName_)[m_scanCount_.get(currentFileName_)] = sx;
    m_scanCount_.put(currentFileName_, m_scanCount_.get(currentFileName_)+1);
  }

  public void AddHeader(CgScanHeader hx) throws CgException
  {
    currentFileName_ = NO_NAME_PLACEHOLDER;
    if (hx.fileName!=null) currentFileName_ = hx.fileName; 
    headerHash_.put(currentFileName_, hx);
    scanHash_.put(currentFileName_, new CgScan[hx.ScanCount]);
    m_scanCount_.put(currentFileName_, 0);
  }

  public void setStartStopHeader(CgScanHeader hx) throws CgException
  {
    if (!headerHash_.containsKey(currentFileName_)) throw new CgException("No header defined");
    headerHash_.get(currentFileName_).StartTime=hx.StartTime;
    headerHash_.get(currentFileName_).EndTime = hx.EndTime;
  }

  public CgScan getLastBaseScan()
  {
    for (int i=(m_scanCount_.get(currentFileName_)-1); i>-1; i--){
      if (scanHash_.containsKey(currentFileName_) && scanHash_.get(currentFileName_)[i].MsLevel==1)
        return scanHash_.get(currentFileName_)[i];
    }
    return null;
  }

  public void addParentFileName(String fileName) throws CgException
  {
    if (this.currentFileName_.equalsIgnoreCase(NO_NAME_PLACEHOLDER)){
      CgScanHeader header = headerHash_.get(currentFileName_);
      header.fileName = fileName;
      headerHash_.remove(currentFileName_);
      CgScan[] scans = scanHash_.get(currentFileName_);
      scanHash_.remove(currentFileName_);
      Integer count = m_scanCount_.get(currentFileName_);
      m_scanCount_.remove(currentFileName_);
      currentFileName_ = fileName;
      header.fileName = fileName;
      headerHash_.put(currentFileName_, header);
      scanHash_.put(currentFileName_, scans);
      m_scanCount_.put(currentFileName_, count);
    }
  }

  public float getLowerThreshold()
  {    
    return this.lowerThresholdTol_;
  }

  public float getUpperThreshold()
  {
    return this.upperThresholdTol_;
  }
  
  /**
   * inits the storage hashes for the read scans
   */
  private void initStoreHashes(){
    m_scanCount_ = new Hashtable<String,Integer>();
    scanHash_ = new Hashtable<String,CgScan[]>();
    headerHash_ = new Hashtable<String,CgScanHeader>();
  }

  /** translates the read scans to the chrom file or a subset chrom file*/
  public void writeToChrom() throws IOException{

    Hashtable<String,Vector<Hashtable<Integer,Float>>> msmsRetentionTimes = buildHeaderInformation();
    Hashtable<String,Vector<LinkedHashMap<String,Vector<MsMsScan>>>> msmsScans = groupTheMsMsScans(msmsRetentionTimes);
    
    Hashtable<String,BufferedOutputStream> stream = new Hashtable<String,BufferedOutputStream>();
    Hashtable<String,DataOutputStream> streamIndex = new Hashtable<String,DataOutputStream>();
    Hashtable<String,float[][]> intensityValuesSection = new Hashtable<String,float[][]>();
    Hashtable<String,Vector<BufferedOutputStream>> msmsChromStreams = new Hashtable<String,Vector<BufferedOutputStream>>();
    Hashtable<String,Vector<DataOutputStream>> msmsIndexStreams = new Hashtable<String,Vector<DataOutputStream>>();
    Hashtable<String,Hashtable<Integer,Integer>> mzIndizes = new Hashtable<String,Hashtable<Integer,Integer>>();
    Hashtable<String,Hashtable<Integer,Long>> bytesIndices = new Hashtable<String,Hashtable<Integer,Long>>();
    
    float floatOfLowerMz = 0f;
    float floatOfUpperMz = 100000f;
    if (lowerThreshInt_!=null && lowerThreshInt_>0){
      floatOfLowerMz = lowerThreshInt_.floatValue()/(float)this.multiplicationFactorForInt_;
      floatOfUpperMz = upperThreshInt_.floatValue()/(float)this.multiplicationFactorForInt_;
    }
    
    for (String key : headerHash_.keySet()){
      stream.put(key, new BufferedOutputStream(new FileOutputStream(dir_+chromFileName_.get(key))));
      streamIndex.put(key,new DataOutputStream(new BufferedOutputStream(new FileOutputStream(dir_+indexFileName_.get(key)))));
      // this is for the MSMS
      msmsChromStreams.put(key, new Vector<BufferedOutputStream>());
      msmsIndexStreams.put(key, new Vector<DataOutputStream>());
      mzIndizes.put(key, new Hashtable<Integer,Integer>());
      bytesIndices.put(key, new Hashtable<Integer,Long>());
      if (msms_&&!msmsInSeveralFiles_){
        for (int i=2; i<=highestMsLevel_;i++){
          BufferedOutputStream streamChrom2 = new BufferedOutputStream(new FileOutputStream(dir_+chromFileName_.get(key)+String.valueOf(i)));
          msmsChromStreams.get(key).add(streamChrom2);
          DataOutputStream streamIndex2 = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(dir_+indexFileName_.get(key)+String.valueOf(i))));
          msmsIndexStreams.get(key).add(streamIndex2);
          int currentResDiffToLowest = 0;
          if (lowerThreshInt_!=null && lowerThreshInt_>0){
            currentResDiffToLowest = (lowerThreshInt_-lowestMz_)/lowestResolution_;
            boolean isFullMultiplicationFactor = true;
            if (currentResDiffToLowest%(elementsForBatchCalculation_*lowestResolution_)!=0) isFullMultiplicationFactor = false;
            currentResDiffToLowest =  (currentResDiffToLowest/(elementsForBatchCalculation_*lowestResolution_))*(elementsForBatchCalculation_*lowestResolution_);
            if (!isFullMultiplicationFactor) currentResDiffToLowest += elementsForBatchCalculation_*lowestResolution_;
          }
          mzIndizes.get(key).put(i-2, currentResDiffToLowest);
          bytesIndices.get(key).put(i-2, 0l);
        }
      }
    }
    writeChrom2(msmsChromStreams,msmsIndexStreams,msmsScans,mzIndizes,bytesIndices,floatOfLowerMz,floatOfUpperMz);
    
    
    Hashtable<String,Long> bytesIndex = new Hashtable<String,Long>();
    for (String key : headerHash_.keySet()){
      bytesIndex.put(key, 0l);
      intensityValuesSection.put(key, createIntensityArray(m_scanCount_.get(key),elementsForBatchCalculation_));
    }

    float intensity;
    if (((this.lowerThreshInt_)/lowestResolution_) % elementsForBatchCalculation_!=0){
      int startValue = this.lowerThreshInt_-((this.lowerThreshInt_-lowestMz_) % (elementsForBatchCalculation_*lowestResolution_));
      for (String key : headerHash_.keySet()){
        for (int j=0; j<m_scanCount_.get(key); j++){
          scanHash_.get(key)[j].fillIntensitiyArray(intensityValuesSection.get(key)[j],(Float.parseFloat((String.valueOf(startValue)))/(float)multiplicationFactorForInt_),(Float.parseFloat((String.valueOf(startValue+elementsForBatchCalculation_*lowestResolution_)))/(float)multiplicationFactorForInt_),elementsForBatchCalculation_);
        }
      }      
    }
    for (int i=this.lowerThreshInt_; i<this.upperThreshInt_; i+=lowestResolution_){
      if (((i-lowestMz_)/lowestResolution_) % elementsForBatchCalculation_==0){
        for (String key : headerHash_.keySet()){
          for (int j=0; j<m_scanCount_.get(key); j++){
            scanHash_.get(key)[j].fillIntensitiyArray(intensityValuesSection.get(key)[j],(Float.parseFloat((String.valueOf(i)))/(float)multiplicationFactorForInt_),(Float.parseFloat((String.valueOf(i+elementsForBatchCalculation_*lowestResolution_)))/(float)multiplicationFactorForInt_),elementsForBatchCalculation_);
          }
        }
      }
      for (String key : headerHash_.keySet()){
        Vector<Integer> scanNumbers = new Vector<Integer>();
        Vector<Float> intensities = new Vector<Float>();
        for (int j=0; j<m_scanCount_.get(key); j++){
          intensity = intensityValuesSection.get(key)[j][((i-lowestMz_)/lowestResolution_) % elementsForBatchCalculation_];
          if (intensity>0){
            scanNumbers.add(new Integer(j));
            intensities.add(intensity);
          }
        }
        String chromString = "";
        if (scanNumbers.size()>0){
          ByteBuffer buffer = ByteBuffer.allocate(scanNumbers.size() * 2 * 4);
          for (int j=0;j!=scanNumbers.size();j++){
            buffer.putInt(scanNumbers.get(j).intValue());
            buffer.putFloat(intensities.get(j).floatValue());
          }
          chromString = String.valueOf(Base64.encode(buffer.array()));
        }
        byte[] bytesToWrite = (chromString+"\n").getBytes();
        stream.get(key).write(bytesToWrite);
        long bytesIdx = bytesIndex.get(key);
        if (((i-lowestMz_)/lowestResolution_) % CgDefines.numberOfEntriesForIndex==0){
          streamIndex.get(key).writeInt((i-lowestMz_)/lowestResolution_);
          streamIndex.get(key).writeLong(bytesIdx);
        }
        bytesIdx+=bytesToWrite.length;
        bytesIndex.put(key, bytesIdx);
      }
    }
    
    numberOfMs2Scans_ = new Hashtable<String,Hashtable<Integer,Integer>>();
    for (String key : headerHash_.keySet()){
      stream.get(key).close();
      stream.get(key).flush();
      streamIndex.get(key).close();
      streamIndex.get(key).flush();
      if (msms_&&!msmsInSeveralFiles_){
        numberOfMs2Scans_.put(key, new Hashtable<Integer,Integer>());
        writeMs2RetentionTimes(dir_,msmsRetentionTimes.get(key),numberOfMs2Scans_.get(key),key);
        for (int i=0; i<highestMsLevel_-1;i++){
          BufferedOutputStream streamChrom2 = msmsChromStreams.get(key).get(i);
          streamChrom2.close();
          streamChrom2.flush();
          DataOutputStream streamIndex2 = msmsIndexStreams.get(key).get(i);
          streamIndex2.close();
          streamIndex2.flush();
        }
      }
    }
    intensityValuesSection = null;
  }
  
  /**
   * 
   * @return the processing hashes for the MSE scans
   */
  private Hashtable<String,Vector<Hashtable<Integer,Float>>> buildHeaderInformation(){
    for (String fileName: headerHash_.keySet()){
      msMsExclusionList_.put(fileName, new Hashtable<Integer,Integer>());
    }
    Hashtable<String,Vector<Hashtable<Integer,Float>>> msmsRetentionTimes = new Hashtable<String,Vector<Hashtable<Integer,Float>>>();
    for (String key : headerHash_.keySet()){
      Vector<Hashtable<Integer,Float>> retsOneMzxml = new Vector<Hashtable<Integer,Float>>();
      for (int i=2; i<=highestMsLevel_;i++) retsOneMzxml.add(new Hashtable<Integer,Float>());
      msmsRetentionTimes.put(key, retsOneMzxml);
    }
    return msmsRetentionTimes;
  }
  
  /**
   * groups the MSn scans correspondingly
   * @param retentionTimes the retention time hash to be filled
   * @return the grouped MSn scans
   */
  @SuppressWarnings("unchecked")
  private Hashtable<String,Vector<LinkedHashMap<String,Vector<MsMsScan>>>> groupTheMsMsScans(Hashtable<String,Vector<Hashtable<Integer,Float>>> retentionTimes){
    Hashtable<String,Vector<LinkedHashMap<String,Vector<MsMsScan>>>> msmsScans = new Hashtable<String,Vector<LinkedHashMap<String,Vector<MsMsScan>>>> ();
    for (String key : this.headerHash_.keySet()){
      Vector<LinkedHashMap<String,Vector<MsMsScan>>> scansOfOneMzXML = new Vector<LinkedHashMap<String,Vector<MsMsScan>>>();
      Hashtable<Integer,Integer> exclusionList = msMsExclusionList_.get(key);
      for (int i=2; i<=highestMsLevel_;i++){
        List<MsMsScan> scansOfOneLevel = new ArrayList<MsMsScan>();
        Hashtable<Integer,Float> rets = retentionTimes.get(key).get(i-2);
        for (CgScan msScan : scanHash_.get(key)){
        //this is because there is too much space allocated in the array
          if (msScan!=null){
            for (CgScan subScan : msScan.getFullSubScans()){
              //the msMsExclusionList is necessary because the scans are read overlapping, if there are several iterations
              if (subScan.MsLevel==i && !exclusionList.containsKey(subScan.Num)){
                scansOfOneLevel.add(((MsMsScan)subScan));
                rets.put(subScan.Num, subScan.RetentionTime);
                exclusionList.put(subScan.Num, subScan.Num);
              }
            }
          }
        }
        //new the scans are sorted according to their precursor mass
        Collections.sort(scansOfOneLevel,new GeneralComparator("at.tugraz.genome.maspectras.quantification.MsMsScan", "getPrecursorMz", "java.lang.Float"));
        LinkedHashMap<String,Vector<MsMsScan>> sortedScans = new LinkedHashMap<String,Vector<MsMsScan>>();
        for (MsMsScan scan : scansOfOneLevel){
          Vector<MsMsScan> samePrecursorMass = new Vector<MsMsScan>();
          if (sortedScans.containsKey(scan.getPrecursorMzAsString())) samePrecursorMass = sortedScans.get(scan.getPrecursorMzAsString());
          samePrecursorMass.add(scan);
          if (samePrecursorMass.size()>1) Collections.sort(samePrecursorMass,new GeneralComparator("at.tugraz.genome.maspectras.quantification.MsMsScan", "getNum", "java.lang.Integer"));
          sortedScans.put(scan.getPrecursorMzAsString(), samePrecursorMass);
        }
        scansOfOneMzXML.add(sortedScans);
      }
      msmsScans.put(key, scansOfOneMzXML);
    }
    return msmsScans;
  }
  
  /** reserves space for an intensity array of a certain size -
   * previously more was done with this array, now the same command may be executed as well directly in the code*/
  private float[][] createIntensityArray(int x, int y){
    float[][] array = new float[x][y];
    return array;
  }

  /**
   * writes the MSn output files (or a subset of them)
   * @param msmsChromStreams the output streams for the chrom2/3 etc. files 
   * @param msmsIndexStreams the output streams for the index files (idx2, idx3, etc.)
   * @param msmsScans the MSn scans
   * @param mzIndices the line index numbers
   * @param bytesIndices the bytes to skip belonging to the line index numbers
   * @param lowerThreshold the lower m/z threshold in integer format
   * @param upperThreshold the upper m/z threshold in integer format
   * @throws IOException thrown if there is a problem with the writing procedure
   */
  private void writeChrom2(Hashtable<String,Vector<BufferedOutputStream>> msmsChromStreams, Hashtable<String,Vector<DataOutputStream>> msmsIndexStreams, Hashtable<String,Vector<LinkedHashMap<String,Vector<MsMsScan>>>> msmsScans,Hashtable<String,Hashtable<Integer,Integer>> mzIndices,Hashtable<String,Hashtable<Integer,Long>> bytesIndices,
      Float lowerThreshold, Float upperThreshold) throws IOException{
    for (String key:headerHash_.keySet()){
      Hashtable<Integer,Integer> mzIndicesOneMzxml  = mzIndices.get(key);
      Hashtable<Integer,Long> byteIndicesOneMzxml = bytesIndices.get(key);
      for (int i=0;i!=msmsChromStreams.get(key).size();i++){
        BufferedOutputStream stream = msmsChromStreams.get(key).get(i);
        DataOutputStream streamIndex = msmsIndexStreams.get(key).get(i);
        LinkedHashMap<String,Vector<MsMsScan>> scans = msmsScans.get(key).get(i);
        int mzIndex = mzIndices.get(key).get(i);
        long bytesIndex = bytesIndices.get(key).get(i);
        for (String scanKey : scans.keySet()){
          Float mzValue = new Float(scanKey);
          if (lowerThreshold!=null && mzValue<lowerThreshold) continue;
        // there are the same amount of index entries as in the normal chrom format
        // e.g. if you have a file from 400-1800Da, a multiplicationFactorForInt_=1000 and a lowestResolution_=2;
        // then you would have (1800-400)*(1000/2)=700 000 lines in the chrom file (because the lowest resolution = 0.002Da)
        // if then CgDefines.numberOfEntriesForIndex = 1000, the index entries must be 700000/1000 = 700;
        // if you would use a lowestResolution_=1, the index file must contain 1400 entries.
        // this calculates the difference in lines to the last lowest mass value (in our case 400Da);
        // e.g. if you have a mass value 402.51Da and the settings with lowestResolution_=2 the calculation would be
        // (402.51*1000-400000)/2=1255 line in the hypothetical chrom file; that means this entry would be after the 2nd index entry (count starts with zero)
          int currentResDiffToLowest = (Math.round(Calculator.roundFloat((Float.parseFloat(scanKey))*(float)this.multiplicationFactorForInt_, 0,BigDecimal.ROUND_UP))-lowestMz_)/lowestResolution_;
        //if the current entry exceeds the value of the next index entry -> write the index entry and increase the index-value
          while (currentResDiffToLowest>=mzIndex && (mzIndex*lowestResolution_)<(upperThreshInt_-lowestMz_)){
            streamIndex.writeInt(mzIndex);
            streamIndex.writeLong(bytesIndex);
            mzIndex+=CgDefines.numberOfEntriesForIndex;
          }
          if (upperThreshold!=null && mzValue>=upperThreshold) continue;
          boolean foundOneScan = false;
          String toWrite = ">"+scanKey+"\n";
          for (MsMsScan scan : scans.get(scanKey)){
            ByteBuffer buffer = ByteBuffer.allocate(scan.Scan.length * 2 * 4);
            for (int j=0; j<scan.Scan.length; j++){
              buffer.putFloat(scan.Scan[j][0]);
              buffer.putFloat(scan.Scan[j][1]);            
            }
            if (scan.Scan.length>0){
              toWrite+=scan.Num+" "+String.valueOf(Base64.encode(buffer.array()))+"\n";
              foundOneScan=true;
            }  
          }
          if (foundOneScan){
            byte[] bytesToWrite = toWrite.getBytes();
            stream.write(bytesToWrite);
            bytesIndex+=bytesToWrite.length;
          }
        }
        //this is for writing the MS2 index in iterations where no MS/MS spectra were found
        if (this.msms_ && scans.size()==0){
          int currentResDiffToLowest = upperThreshInt_-lowestMz_;
        //if the current entry exceeds the value of the next index entry -> write the index entry and increase the index-value
          while (currentResDiffToLowest>=mzIndex && (mzIndex*lowestResolution_)<(upperThreshInt_-lowestMz_)){
            streamIndex.writeInt(mzIndex);
            streamIndex.writeLong(bytesIndex);
            mzIndex+=CgDefines.numberOfEntriesForIndex;
          }

        }
        byteIndicesOneMzxml.put(i, bytesIndex);
        mzIndicesOneMzxml.put(i, mzIndex);
      }
      mzIndices.put(key, mzIndicesOneMzxml);
      bytesIndices.put(key, byteIndicesOneMzxml);
    }
  }

  /**
   * writess the MSn retention time information
   * @param directory the output directory
   * @param msmsRetentionTimes the MSn retention times
   * @param numberOfMs2Scans the respective amount of MSn scans
   * @param fileName key for merged files
   * @throws IOException
   */
  private void writeMs2RetentionTimes(String directory, Vector<Hashtable<Integer,Float>> msmsRetentionTimes, Hashtable<Integer,Integer> numberOfMs2Scans, String fileName)throws IOException{
    for (int i=2; i<=highestMsLevel_;i++){
      Hashtable<Integer,Float> rets = msmsRetentionTimes.get(i-2);
      DataOutputStream streamRetTime = null;
      List<Integer> scans = new ArrayList<Integer>();
      Set<Integer> keys = rets.keySet();
      for (Integer key : keys)scans.add(key);
      Collections.sort(scans);
      try{
        int count = 0;
        streamRetTime = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(directory+this.retentionTimeFileName_.get(fileName)+String.valueOf(i))));
        for (Integer key : scans){
          streamRetTime.writeInt(key);
          streamRetTime.writeFloat(rets.get(key));
          count++;
        }
        numberOfMs2Scans.put(i, count);
      }catch (IOException iox){
        iox.printStackTrace();
      }finally{
        if (streamRetTime!=null){
          streamRetTime.close();
          streamRetTime.flush();
        }
      }
    }
  }
  
  /**
   * sets general properties for the chrom translation
   * @param elementsForBatchCalculation the arraysize for reading m/z ranges in batch mode
   * @param lowestMz the lowest possible m/z value
   * @param highestMz the highest possible m/z value
   * @param chromFileName the names of the chrom files
   * @param indexFileName the names of the index files
   * @param retentionTimeFileName the names of the retention time files
   */
  public void setStaticInformation(int elementsForBatchCalculation, int lowestMz, int highestMz, 
      Hashtable<String,String> chromFileName, Hashtable<String,String> indexFileName, Hashtable<String,String> retentionTimeFileName)
  {
    this.elementsForBatchCalculation_ = elementsForBatchCalculation;
    this.lowestMz_ = lowestMz;
    this.highestMz_ = highestMz;
    for (String key : chromFileName.keySet()){
      String name = StringUtils.getJustFileName(chromFileName.get(key));
      this.chromFileName_.put(key, name);
      name = StringUtils.getJustFileName(indexFileName.get(key));
      this.indexFileName_.put(key, name);
      name = StringUtils.getJustFileName(retentionTimeFileName.get(key));
      this.retentionTimeFileName_.put(key, name);
    }    
  }
  
  /**
   * in the case of one iteration, and only one thread - the file is already read in the readHeaderInformation
   * this method passes these information to this object 
   * @param headerHash the header information
   * @param scanHash the individual scans
   * @param m_scanCount actual Scan Count (used as index)
   */
  public void setReadXmlContent(Hashtable<String,CgScanHeader> headerHash, Hashtable<String,CgScan[]> scanHash,Hashtable<String,Integer> m_scanCount){
    this.headerHash_ = headerHash;
    this.scanHash_ = scanHash;
    this.m_scanCount_ = m_scanCount;
  }
  
  public Hashtable<String,Integer> getm_scanCount()
  {
    return m_scanCount_;
  }
  
  /**
   * writes the retention-time file
   * @throws IOException
   */
  public void writeRetentionTimeFile(Hashtable<String,String> retentionTimeFileName) throws IOException{
    for (String key : retentionTimeFileName.keySet()){
      DataOutputStream streamRetTime = null;
      try{
        streamRetTime = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(retentionTimeFileName.get(key))));
        for (int i=0; i<m_scanCount_.get(key); i++){
          streamRetTime.writeInt(this.scanHash_.get(key)[i].Num);
          streamRetTime.writeFloat(this.scanHash_.get(key)[i].RetentionTime);
        }
      }catch (IOException iox){
        iox.printStackTrace();
      }finally{
        if (streamRetTime!=null){
          streamRetTime.close();
          streamRetTime.flush();
        }
      }
    }
  }

  /** returns the amount of MS2 scans*/
  public Hashtable<String,Hashtable<Integer,Integer>> getNumberOfMs2Scans()
  {
    return numberOfMs2Scans_;
  }
  
  /**
   * 
   * @return true when the thread terminated
   */
  public boolean finished(){
    return this.finished_;
  }

  /**
   * 
   * @return the error message in case of an error
   */
  public String getErrorString(){
    return this.errorString_;
  }
  
  /**
   *  cleans the storage hashes to release memory space
   */
  @SuppressWarnings("unused")
  public void cleanUp(){
    for (String key : this.headerHash_.keySet()){
      CgScan[] scans = scanHash_.get(key);
      for (CgScan scan : scans){
        scan = null;
      }
      CgScanHeader header = headerHash_.get(key);
      header = null;
      Hashtable<Integer,Integer> exclusion = msMsExclusionList_.get(key);
      exclusion = null;
      Hashtable<Integer,Integer> numberOfMs2 = numberOfMs2Scans_.get(key);
      numberOfMs2 = null;
    }
    scanHash_ = null;
    headerHash_ = null;
    msMsExclusionList_ = null;
    numberOfMs2Scans_ = null;    
  }
  
  /**
   * calculates the integer format of an m/z value
   * @param mzValue the m/z value to be converted in integer
   * @param multiplicationFactorForInt
   * @return the multiplication factor to be used, to create integer values out of the m/z float values
   */
  private static int getCorrespondingIntValue(float mzValue, int multiplicationFactorForInt){
    return (int)Calculator.roundFloat((mzValue*(float)multiplicationFactorForInt),0,BigDecimal.ROUND_DOWN);
  }

  
}
