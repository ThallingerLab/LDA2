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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.Timer;
import java.util.TimerTask;

import at.tugraz.genome.lda.utils.RangeInteger;
import at.tugraz.genome.maspectras.GlobalConstants;
import at.tugraz.genome.maspectras.quantification.CgDefines;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgScan;
import at.tugraz.genome.maspectras.quantification.CgScanHeader;
import at.tugraz.genome.maspectras.quantification.ChromatogramReader;
import at.tugraz.genome.maspectras.utils.StringUtils;
import at.tugraz.genome.util.BioUtilsConstants;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class RawToChromTranslator implements AddScan
{
  
  /** actual Scan Count (used as index) */
  private Hashtable<String,Integer> m_scanCount;
  /** actual polarity scan count (used as index) */
  private Hashtable<String,Hashtable<Integer,Integer>> polarity_scanCount_;
  /** The mzXml File Name */
  private String m_fileName;
  /** the lowest possible mz value*/
  private int lowestMz;
  /** the highest possible mz value*/
  private int highestMz;
  /** the name of the chromatography file(s)*/
  private Hashtable<String,String[]> chromFileName_;
  /** the name of the index file(s)*/
  private Hashtable<String,String[]> indexFileName_;
  /** the name of the header file(s)*/
  private Hashtable<String,String[]> headerFileName_;
  /** the name of the retention-time file(s)*/
  private Hashtable<String,String[]> retentionTimeFileName_;
  
  /** the lower and upper m/z thresholds for every thread*/
  private Hashtable<Integer,Hashtable<Integer,RangeInteger>> threadBoundaries_;
  
  /** the reader*/
  protected XmlSpectraReader m_reader;
  /** the available threads for chrom writing*/
  private Hashtable<Integer,RawToChromThread> translators_;
  
  /** the multiplication factor to be used, to create integer values out of the m/z float values*/
  private int multiplicationFactorForInt_;
  /** the lowest resolution distance value by a given multiplicationFactorForInt_*/ 
  private int lowestResolution_;
  /** the arraysize for reading m/z ranges in batch mode*/
  private int elementsForBatchCalculation_;
  
  /** the type of the data file, currently, only mzXML is supported*/
  private String fileType_;
  /** the highest number of MBytes that can be translated in one chrom translation iteration*/
  private int maxMBForChromTranslation_;
  /** the amount of available threads*/
  private int numberOfThreads_;
  /** allow reading of MSn information*/
  private boolean msms_;
  /** the highest MSn level in the data file*/
  private int highestMsLevel_;
  /** contains the amount of MSn scans - first key: file key; second key: polarity; third key: MS-level*/
  private Hashtable<String,Hashtable<Integer,Hashtable<Integer,Integer>>> numberOfMs2Scans_;
  
  /** is MSn information distributed in several files, such as in Waters MSE data*/
  private boolean msmsInSeveralFiles_;
  /** the highest MSn level, if MSn information is distributed in several files, such as in Waters MSE data*/
  private int msmsSeveralFilesHighestLevel_;
  
  /** the lowest m/z values of the various MSn levels in integer format*/
  private Hashtable<Integer,Integer> lowestMzs_;
  /** the highest m/z values of the various MSn levels in integer format*/
  private Hashtable<Integer,Integer> highestMzs_;
  /** the number of scans for each MSn level*/
  private Hashtable<Integer,Hashtable<String,Integer>> scanCounts_;
  /** the number of scans for each MSn level, classified in polarities*/
  private Hashtable<Integer,Hashtable<String,Hashtable<Integer,Integer>>> polarityScanCounts_;
  
  /** the header information*/
  private Hashtable<String,CgScanHeader> headerHash_;
  /** the individual scans*/
  private Hashtable<String,CgScan[]> scanHash_;
  /** the highest MSn level depending on separate data files in one physical file, and the polarity*/
  private Hashtable<String,Hashtable<Integer,Integer>> highestMsLevelHash_;
  
  /** the data file name (in case of merged data) that is currently read by the reader*/ 
  private String currentFileName_;
  
  /** if there are no merged MS runs, this value is used as key for the hashtables*/ 
  private final static String NO_NAME_PLACEHOLDER = "";
  
  /** the timer for multi-threaded writing of chrom files*/
  private Timer timer_;
  /** if all threads are finished, one chrom iteration is finished*/
  private boolean oneIterationFinished_ = false;
  /** the error string, if any thread encounters an error*/
  private String errorString_ = null;
  /** the current status of the single threads*/
  private Hashtable<Integer,Integer> quantStatus_;
  
  /** thread status: waiting for job*/
  private final static int STATUS_WAITING = 0;
  /** thread status: excecuting the job*/
  private final static int STATUS_CALCULATING = 1;
  /** thread status: thread is finished*/
  private final static int STATUS_FINISHED = 2;
  
  /** the buffer size for appending files*/
  private static final int BUFFER_SIZE = 1024 * 4;

  /** the lower m/z threshold for data reading*/
  private float lowerThreshold_ = 0;
  /** the upper m/z threshold for data reading*/
  private float upperThreshold_ = 1000000;
  /** the required chrom iterations, due to memory consumption*/
  private int numberOfIterations_;
  
  /** was polarity switching used*/
  private boolean polaritySwitching_ = false;
  
  /**
   * constructor for the translation
   * @param mzXmlPath the original file
   * @param fileType the type of the file
   * @param maxMBForChromTranslation the highest number of MBytes that can be translated in one chrom translation iteration
   * @param multiplicationFactorForInt the multiplication factor to be used, to create integer values out of the m/z float values
   * @param lowestResulution the multiplication factor to be used, to create integer values out of the m/z float values
   */
  public RawToChromTranslator(String mzXmlPath,String fileType, int maxMBForChromTranslation, int multiplicationFactorForInt,
      int lowestResulution)
  {
   this(mzXmlPath,fileType,maxMBForChromTranslation);
   multiplicationFactorForInt_ = multiplicationFactorForInt;
   lowestResolution_ = lowestResulution;
   elementsForBatchCalculation_ = multiplicationFactorForInt_/lowestResolution_;
  }
  
  /**
   * constructor for the translation
   * @param mzXmlPath the original file
   * @param fileType the type of the file
   * @param maxMBForChromTranslation the highest number of MBytes that can be translated in one chrom translation iteration
   */
  public RawToChromTranslator(String mzXmlPath,String fileType, int maxMBForChromTranslation)
  {
    multiplicationFactorForInt_ = CgDefines.mzMultiplicationFactorForInt;
    lowestResolution_ = CgDefines.lowestResolution;
    elementsForBatchCalculation_ = CgDefines.elementsForBatchCalculation;
    m_fileName = mzXmlPath;
    fileType_ = fileType;
    maxMBForChromTranslation_ = maxMBForChromTranslation;
    msms_ = false;
    highestMsLevelHash_ = new Hashtable<String,Hashtable<Integer,Integer>>();    
    initStoreHashes();
    highestMsLevel_ = 1;
    numberOfThreads_ = 1;
  }

  /**
   * constructor for the translation
   * @param mzXmlPath the original file
   * @param fileType the type of the file
   * @param maxMBForChromTranslation the highest number of MBytes that can be translated in one chrom translation iteration
   * @param numberOfThreads the number of threads that shall be used for the translation
   * @param multiplicationFactorForInt the multiplication factor to be used, to create integer values out of the m/z float values
   * @param lowestResulution the multiplication factor to be used, to create integer values out of the m/z float values
   * @param msms allow reading of MSn information
   */
  public RawToChromTranslator(String mzXmlPath,String fileType, int maxMBForChromTranslation, int numberOfThreads,
      int multiplicationFactorForInt, int lowestResulution, boolean msms)
  {
   this(mzXmlPath,fileType,maxMBForChromTranslation,multiplicationFactorForInt,lowestResulution);
   msms_ = msms;
   numberOfThreads_ = numberOfThreads;
  }
  
  public void AddHeader(CgScanHeader hx) throws CgException
  {
    currentFileName_ = NO_NAME_PLACEHOLDER;
    if (hx.fileName!=null) currentFileName_ = hx.fileName; 
    headerHash_.put(currentFileName_, hx);
    scanHash_.put(currentFileName_, new CgScan[hx.ScanCount]);
    m_scanCount.put(currentFileName_, 0);
    
    Hashtable<Integer,Integer> polarityScans = new Hashtable<Integer,Integer>();
    polarityScans.put(CgDefines.POLARITY_NO, 0);
    polarityScans.put(CgDefines.POLARITY_POSITIVE, 0);
    polarityScans.put(CgDefines.POLARITY_NEGATIVE, 0);
    polarity_scanCount_.put(currentFileName_, polarityScans);    
  }

  public void AddScan(CgScan sx) throws CgException
  {
    if (!scanHash_.containsKey(currentFileName_)) throw new CgException("m_scans Array not allocated");
    if (scanHash_.get(currentFileName_).length==m_scanCount.get(currentFileName_)) return;
    scanHash_.get(currentFileName_)[m_scanCount.get(currentFileName_)] = sx;
    m_scanCount.put(currentFileName_, m_scanCount.get(currentFileName_)+1);
    polarity_scanCount_.get(currentFileName_).put(sx.getPolarity(),polarity_scanCount_.get(currentFileName_).get(sx.getPolarity())+1);
  }

  public void setStartStopHeader(CgScanHeader hx) throws CgException
  {
    if (!headerHash_.containsKey(currentFileName_)) throw new CgException("No header defined");
    headerHash_.get(currentFileName_).StartTime=hx.StartTime;
    headerHash_.get(currentFileName_).EndTime = hx.EndTime;
  }

  public CgScan getLastBaseScan(){
    for (int i=(m_scanCount.get(currentFileName_)-1); i>-1; i--){
      if (scanHash_.containsKey(currentFileName_) && scanHash_.get(currentFileName_)[i].MsLevel==1)
        return scanHash_.get(currentFileName_)[i];
    }
    return null;
  }
  
  public void addParentFileName(String fileName) throws CgException{
    if (this.currentFileName_.equalsIgnoreCase(NO_NAME_PLACEHOLDER)){
      CgScanHeader header = headerHash_.get(currentFileName_);
      header.fileName = fileName;
      headerHash_.remove(currentFileName_);
      CgScan[] scans = scanHash_.get(currentFileName_);
      scanHash_.remove(currentFileName_);
      Integer count = m_scanCount.get(currentFileName_);
      m_scanCount.remove(currentFileName_);
      Hashtable<Integer,Integer> polarityCount = polarity_scanCount_.get(currentFileName_);
      polarity_scanCount_.remove(currentFileName_);
      currentFileName_ = fileName;
      header.fileName = fileName;
      headerHash_.put(currentFileName_, header);
      scanHash_.put(currentFileName_, scans);
      m_scanCount.put(currentFileName_, count);
      polarity_scanCount_.put(currentFileName_, polarityCount);
    }
  }
  
  /**
   * inits the chrom translation process
   * @throws CgException the exception if anything is wrong
   */
  public void translateToChromatograms() throws CgException{
    msmsInSeveralFiles_ = false;
    msmsSeveralFilesHighestLevel_ = 1;
    lowestMzs_ = new Hashtable<Integer,Integer>();
    highestMzs_  = new Hashtable<Integer,Integer>();
    scanCounts_ = new Hashtable<Integer,Hashtable<String,Integer>>();
    polarityScanCounts_ = new Hashtable<Integer,Hashtable<String,Hashtable<Integer,Integer>>>();
    oneIterationFinished_ = false;
    errorString_ = null;
    if (msms_){
      int msLevel = 2;
      while ((new File(m_fileName+String.valueOf(msLevel))).exists()){
        msmsInSeveralFiles_ = true;
        msmsSeveralFilesHighestLevel_ = msLevel;
        msLevel++;
      }
      for (int i=msmsSeveralFilesHighestLevel_;i!=0;i--)translateToChromatograms(i); 
    }else{
      translateToChromatograms(1);
    }
    try {
      this.writeHeaderFile();
    }
    catch (IOException e) {
      e.printStackTrace();
      throw new CgException(e.getMessage());
    }

  }

  /**
   * 
   * @return the original file names (in case of merged files)
   */
  public Set<String> getOriginalFileNames(){
    return this.headerHash_.keySet();
  }

  /**
   * starts the translation process
   * @param msLevel for MSE files, the file is stored in several levels (files)
   * @throws CgException the exception if anything is wrong
   */
  private void translateToChromatograms(int msLevel) throws CgException{
    long time = System.currentTimeMillis();
    this.readHeaderInformation(msLevel);
        
    if (this.numberOfIterations_>1 || this.numberOfThreads_>1){
      m_scanCount = new Hashtable<String,Integer>();
    }
    String suffix = "";
    if (msLevel>1) suffix  = String.valueOf(msLevel);
    this.initTranslatorObjects();
    for (int i=0; i!=this.numberOfIterations_; i++){
      System.out.println("Starting iteration: "+(i+1));
      this.quantStatus_ = new Hashtable<Integer,Integer>();
      for (int j=0; j!=this.numberOfThreads_; j++){
        RangeInteger threshold = this.getLowerUpperThreshold(i, j);
        String[] directories = new String[2];
        if (polaritySwitching_){
          directories[0] = getCorrespondingThreadDirectory(CgDefines.POLARITY_POSITIVE, i, j);
          directories[1] = getCorrespondingThreadDirectory(CgDefines.POLARITY_NEGATIVE, i, j);
        }else{
          directories[0] = getCorrespondingThreadDirectory(CgDefines.POLARITY_NO, i, j);
          directories[1] = null;
        }
        RawToChromThread thread = new RawToChromThread(this.translators_.get(j));
        thread.setRequiredInformation(directories,threshold.getStart(),threshold.getStop()); 
        this.translators_.put(j, thread);
        this.quantStatus_.put(j, STATUS_WAITING);
      }
      AddScan[] adders = getAdders();
      m_reader.setAdders(adders);
      if (this.numberOfIterations_>1 || this.numberOfThreads_>1){
        this.m_reader.ReadFile(m_fileName+suffix);
      }else{
        this.translators_.get(0).setReadXmlContent(headerHash_,scanHash_,m_scanCount,polarity_scanCount_);
      }
    
      if (this.numberOfThreads_>1){
        timer_ = new java.util.Timer();
        ThreadSupervisor supervisor = new ThreadSupervisor(msLevel,i);
        timer_.schedule(supervisor, 10, 100);
        while (!oneIterationFinished_){
          try {
            this.wait(1000);
          }
          catch (Exception ex) {
          }
        }
        try {
          this.wait(500);
        }
        catch (Exception ex) {
        }
        oneIterationFinished_ = false;
        this.timer_.cancel();
        this.timer_ = null;
        supervisor = null;
      }else{
        try{
          RawToChromThread singleThread = translators_.get(0);
          String dir = this.m_fileName.substring(0,this.m_fileName.lastIndexOf("."))+".chrom/";
          if (i>0){
            dir +=String.valueOf(i)+"/";
            File dirFile = new File(dir);
            dirFile.mkdir();
          }
          singleThread.writeToChrom();   
          if (i==0){
            m_scanCount = singleThread.getm_scanCount();
            polarity_scanCount_ = singleThread.getPolarityScanCount();
            scanCounts_.put(msLevel, m_scanCount);
            polarityScanCounts_.put(msLevel, polarity_scanCount_);
            try{
              singleThread.writeRetentionTimeFile(retentionTimeFileName_,polaritySwitching_);
            } catch (IOException iox){
              iox.printStackTrace();
              errorString_ = iox.getMessage();
            }
            this.numberOfMs2Scans_ = singleThread.getNumberOfMs2Scans();
          }
          singleThread.cleanUp();
          System.gc();
        }catch (Exception ex){
          ex.printStackTrace();
          errorString_ = ex.getMessage();
        }
      }
    }
    for (@SuppressWarnings("unused") RawToChromThread  translator : translators_.values()){
      translator = null;
    }
    translators_ = null;

    if (errorString_!=null) throw new CgException(errorString_);
      
    mergeResults();
    try {
      checkHighestMsLevel();
    }
    catch (IOException e) {
      e.printStackTrace();
    }
    System.out.println("Total time: "+((System.currentTimeMillis()-time)/1000)+" secs");
  }

  /**
   * merges the different chrom translations in a single file - when the file was translated in several slices
   */
  private void mergeResults(){
    if (this.numberOfIterations_>1 || this.numberOfThreads_>1){
      numberOfMs2Scans_ = new Hashtable<String,Hashtable<Integer,Hashtable<Integer,Integer>>>();
      mergeResults(0);
      if (polaritySwitching_)
        mergeResults(1);
    }
  }
  
  /**
   * merges the different chrom translations in a single file - when the file was translated in several slices
   * @param filePosition when polarity switching is used, two chrom directories are generated - the file position is which one shall be merged
   */
  private void mergeResults(int filePosition){
    int totalSlices = this.numberOfIterations_*this.numberOfThreads_;    
    try{
      for (String key:headerHash_.keySet()){
        String baseDir = getBaseDir(filePosition);
        if (msms_&&!msmsInSeveralFiles_){
          Hashtable<Integer,Hashtable<Integer,Integer>> msmsScans = new Hashtable<Integer,Hashtable<Integer,Integer>>();
          if (numberOfMs2Scans_.containsKey(key))
            msmsScans = numberOfMs2Scans_.get(key);
          msmsScans.put(filePosition, new Hashtable<Integer,Integer>());
          numberOfMs2Scans_.put(key, msmsScans);
        }
        for (int i=1; i<=highestMsLevel_;i++){
          String levelSuffix = "";
          if (i>1) levelSuffix = String.valueOf(i);
          //merge index
          Hashtable<Integer,Integer> lines = new Hashtable<Integer,Integer>();
          Hashtable<Integer,Long> indices = new Hashtable<Integer,Long>();
          int count = 0;
          long previousLength = 0;
          DataInputStream indexStream = new DataInputStream(new FileInputStream(indexFileName_.get(key)[filePosition]+levelSuffix));
          count = readIndexFile(indexStream,lines,indices,count,previousLength);
          indexStream.close();
          previousLength += (new File (chromFileName_.get(key)[filePosition]+levelSuffix)).length();
          for (int j=1; j!=totalSlices; j++){
            String dir = baseDir+String.valueOf(j)+"/";
            String indexFileName = dir+StringUtils.getJustFileName(indexFileName_.get(key)[filePosition])+levelSuffix;
            DataInputStream in = new DataInputStream(new FileInputStream(indexFileName));
            count = readIndexFile(in,lines,indices,count,previousLength);
            in.close();
            previousLength += (new File(dir+StringUtils.getJustFileName(chromFileName_.get(key)[filePosition])+levelSuffix)).length();
          }
          DataOutputStream streamIndex2 = new DataOutputStream(new FileOutputStream(indexFileName_.get(key)[filePosition]+levelSuffix));
          for (int j=0; j!=count; j++){
            streamIndex2.writeInt(lines.get(j));
            streamIndex2.writeLong(indices.get(j));
          }
          streamIndex2.close();
            //merge retention times
          if (i>1){
            Hashtable<Integer,Float> retentionTimes = new Hashtable<Integer,Float>();
            DataInputStream rttStream = new DataInputStream(new FileInputStream(retentionTimeFileName_.get(key)[filePosition]+levelSuffix));
            readRetentionTimeFile(rttStream,retentionTimes);
            rttStream.close();
            for (int j=1; j!=totalSlices; j++){
              String dir = baseDir+String.valueOf(j)+"/";
              String rttFileName = dir+StringUtils.getJustFileName(retentionTimeFileName_.get(key)[filePosition])+levelSuffix;
              DataInputStream in = new DataInputStream(new FileInputStream(rttFileName));
              readRetentionTimeFile(in,retentionTimes);
              in.close();
            }
            List<Integer> scanNumbers = new ArrayList<Integer>(retentionTimes.keySet());
            Collections.sort(scanNumbers);
            DataOutputStream streamRtt2 = new DataOutputStream(new FileOutputStream(retentionTimeFileName_.get(key)[filePosition]+levelSuffix));
            for (Integer scanNumber : scanNumbers){
              streamRtt2.writeInt(scanNumber);
              streamRtt2.writeFloat(retentionTimes.get(scanNumber));
            }
            streamRtt2.close();
            numberOfMs2Scans_.get(key).get(filePosition).put(i, scanNumbers.size());
          }
          //merge chrom files
          BufferedOutputStream streamChrom2 = new BufferedOutputStream(new FileOutputStream(chromFileName_.get(key)[filePosition]+levelSuffix,true));
          for (int j=1; j!=totalSlices; j++){
            String dir = baseDir+String.valueOf(j)+"/";
            BufferedInputStream in = new BufferedInputStream(new FileInputStream(dir+StringUtils.getJustFileName(chromFileName_.get(key)[filePosition])+levelSuffix));
            RawToChromTranslator.append(in,streamChrom2);
            in.close();
          }
          streamChrom2.close();

        }
        for (int j=1; j!=totalSlices; j++){
          String dir = baseDir+String.valueOf(j)+"/";
          File dirFile = new File(dir);
          for (File file : dirFile.listFiles()){
            file.delete();
          }
          dirFile.delete();
        }
      }
    }catch(IOException iox){
      iox.printStackTrace();
    }    
  }
  
  /**
   * checks if for this chrom file really all MS-levels are present - otherwise, the additional (empty) files are deleted 
   * @throws IOException thrown when there is something wrong with the file access
   */
  private void checkHighestMsLevel() throws IOException{
    for (String key:headerHash_.keySet()){
      Hashtable<Integer,Integer> polarityMsLevels =  highestMsLevelHash_.get(key);
      for (Integer filePosition : polarityMsLevels.keySet()){
        int highestMsLevel = polarityMsLevels.get(filePosition);
        while (highestMsLevel>1){
          String levelSuffix = String.valueOf(highestMsLevel);
          File rttFile = new File(retentionTimeFileName_.get(key)[filePosition]+levelSuffix);
          //the file at this level is not empty, stop the procedure
          if (rttFile.length()!=0l)
            break;
          File indexFile = new File(indexFileName_.get(key)[filePosition]+levelSuffix);
          File chromFile = new File(chromFileName_.get(key)[filePosition]+levelSuffix);
          chromFile.delete();
          indexFile.delete();
          rttFile.delete();
          highestMsLevel--;
        }
        polarityMsLevels.put(filePosition,highestMsLevel);
      }
    }
  }

  /**
   * reads the header information of the corresponding file
   * @param msLevel for MSE files, the file is stored in several levels (files)
   */
  private void readHeaderInformation(int msLevel){
    String suffix = "";
    chromFileName_ = new Hashtable<String,String[]>();
    indexFileName_ = new Hashtable<String,String[]>();
    headerFileName_ = new Hashtable<String,String[]>();
    retentionTimeFileName_ = new Hashtable<String,String[]>();
    if (msLevel>1) suffix = String.valueOf(msLevel);
    m_scanCount = new Hashtable<String,Integer>();
    threadBoundaries_ = new Hashtable<Integer,Hashtable<Integer,RangeInteger>>();
    
    File fileInfo = new File (m_fileName+suffix);
    long fileSize = fileInfo.length();
    numberOfIterations_ = Integer.parseInt(Long.toString(fileSize/(((long)this.maxMBForChromTranslation_)*1024l*1024l)))+1;
    
    System.out.println("Number of iterations: "+numberOfIterations_);
     
    try
    {
      m_reader = null;
      String justFileName = "";
      //check if several mzXML files exist (several MS/MS traces)
      if (fileType_.equalsIgnoreCase("mzXML")){
        AddScan[] adders = new AddScan[1];
        adders[0] = this;
        m_reader = new MzXmlReader(adders, msms_&&!msmsInSeveralFiles_,multiplicationFactorForInt_);
        if (m_fileName.endsWith(".mzXML")) justFileName = m_fileName.substring(0, m_fileName.length()-(".mzXML").length());
        else justFileName = m_fileName; 
      }
/*    these are not updated to the parallelized reading
      if (fileType_.equalsIgnoreCase("mzData")){
        System.out.println("Taking a mzData File");
        m_reader = new CgMzDataReader(this,multiplicationFactorForInt_);
        if (m_fileName.endsWith(".mzData")) justFileName = m_fileName.substring(0, m_fileName.length()-(".mzData").length());
        else justFileName = m_fileName;
      }
      if (fileType_.equalsIgnoreCase("RAW")){
        System.out.println("Taking a XCalibur RAW File");
        m_reader = new CgXCaliburRawReader(this,multiplicationFactorForInt_);
        if (m_fileName.endsWith(".RAW")) justFileName = m_fileName.substring(0, m_fileName.length()-(".RAW").length());
        else justFileName = m_fileName;
      }*/
      if (m_reader!=null){
        if (this.numberOfIterations_>1 || this.numberOfThreads_>1){
          m_reader.ReadFile(m_fileName+suffix,true);
        }else{
          m_reader.ReadFile(m_fileName+suffix);
        }
        polaritySwitching_ = m_reader.usesPolaritySwitching();
        if (this.headerHash_.size()==1){
          String justTheName = StringUtils.getJustFileName(justFileName);
          generateChromDirectoryAndFiles(justFileName,justTheName,headerHash_.keySet().iterator().next(), suffix);
        } else if (this.headerHash_.size()>1){
          String directory = StringUtils.extractDirName(m_fileName);
          for (String name : headerHash_.keySet()){
            String chromDirName = directory+File.separator+name;
            generateChromDirectoryAndFiles(chromDirName,name,name,suffix);
          }
        }
        this.highestMz = m_reader.getHighestMz();
        this.lowestMz = m_reader.getLowestMz();
        lowestMzs_.put(msLevel, lowestMz);
        highestMzs_.put(msLevel, highestMz);
        int nrOfLLowResItems = this.highestMz-this.lowestMz;
        int previousStartLine = 0;
        for (int i=0; i!=numberOfIterations_;i++){
          int lineNumber = ((i+1)*nrOfLLowResItems)/numberOfIterations_;
          if (i==numberOfIterations_-1) lineNumber=this.highestMz-this.lowestMz;
          //this is to generate exact borders for lowestResolutions>1
          else if (lineNumber%this.lowestResolution_!=0){
            if (lineNumber%this.lowestResolution_<(this.lowestResolution_-(lineNumber%this.lowestResolution_)))
              lineNumber-=(lineNumber%this.lowestResolution_);
            else
              lineNumber+=this.lowestResolution_-(lineNumber%this.lowestResolution_);

          }
          int previousThreadLine = previousStartLine;
          Hashtable<Integer,RangeInteger> rangeBoundaries = new  Hashtable<Integer,RangeInteger>();
          for (int j=0; j!=numberOfThreads_; j++){
            int threadLine = ((j+1)*(lineNumber-previousStartLine))/numberOfThreads_+previousStartLine;
            if (j==numberOfThreads_-1) threadLine=lineNumber;
            //this is to generate exact borders for lowestResolutions>1
            else if (threadLine%this.lowestResolution_!=0){
              if (threadLine%this.lowestResolution_<(this.lowestResolution_-(threadLine%this.lowestResolution_)))
                threadLine-=(threadLine%this.lowestResolution_);
              else
                threadLine+=this.lowestResolution_-(threadLine%this.lowestResolution_);
            }
            rangeBoundaries.put(j, new RangeInteger(previousThreadLine+this.lowestMz,new Integer(threadLine+this.lowestMz)));
            previousThreadLine = threadLine;
          }
          threadBoundaries_.put(i, rangeBoundaries);
          previousStartLine = lineNumber;
        }
      }else throw new CgException("No readable file format.");
    }
    catch(Exception ex)
    {
      ex.printStackTrace();
    }
    for (String fileName: headerHash_.keySet()){
      if (headerHash_.get(fileName).highestMSLevel>highestMsLevel_)
        highestMsLevel_ = headerHash_.get(fileName).highestMSLevel;
      Hashtable<Integer,Integer> polarityMsLevels = new Hashtable<Integer,Integer>();
      polarityMsLevels.put(0, headerHash_.get(fileName).highestMSLevel);
      if (polaritySwitching_)
        polarityMsLevels.put(1, headerHash_.get(fileName).highestMSLevel);
      highestMsLevelHash_.put(fileName, polarityMsLevels);
    }
  }
  
  /**
   * inits the thread objects for chrom translation
   */
  private void initTranslatorObjects(){
    translators_ = new Hashtable<Integer,RawToChromThread>();
    for (int i=0; i!=this.numberOfThreads_; i++){
      RawToChromThread thread = new RawToChromThread(this.msms_,this.multiplicationFactorForInt_,this.lowestResolution_,this.highestMsLevel_);
      thread.setStaticInformation(elementsForBatchCalculation_,lowestMz,highestMz,chromFileName_,indexFileName_,
          retentionTimeFileName_);
      this.translators_.put(i,thread);
    }
  }
  
  /**
   * writes the header file
   */
  private void writeHeaderFile() throws IOException {
    for (String key : this.headerFileName_.keySet()){
      if (polaritySwitching_){
        this.writeHeaderFile(headerFileName_.get(key)[0],key,CgDefines.POLARITY_POSITIVE);
        this.writeHeaderFile(headerFileName_.get(key)[1],key,CgDefines.POLARITY_NEGATIVE);
      }else{
        this.writeHeaderFile(headerFileName_.get(key)[0],key,CgDefines.POLARITY_NO);      
      }
    }
//    for (String key : this.headerFileName_.keySet()){
//    }
  }
  
  /**
   * writes the header file - in the case of polarity switching, more than one header file has to be written
   * @param fileName the path to the header file
   * @param key the file key
   * @param polarity the current polarity to be written (as defined in CgDefines)
   * @throws IOException when there is something wrong with the file access
   */
  private void writeHeaderFile(String fileName, String key, int polarity) throws IOException{
    int filePosition = 0;
    if (polarity == CgDefines.POLARITY_NEGATIVE) filePosition = 1;
    Properties props = new Properties();
    props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_MZ_MULTIPLICATION_FACTOR, String.valueOf(multiplicationFactorForInt_));
    props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_LOWEST_RESOLUTION,String.valueOf(lowestResolution_));
    props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_LOWEST_MZ, String.valueOf(this.lowestMz));
    props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_HIGHEST_MZ, String.valueOf(this.highestMz));
    props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_NUMBER_SCANS, String.valueOf(this.getCorrectScanCount(key, polarity)));   
    props.put(BioUtilsConstants.INDEX_HEADER_FILE_NUMBER_OF_ENTRIES_FOR_INDEX, String.valueOf(CgDefines.numberOfEntriesForIndex));
    props.put(BioUtilsConstants.INDEX_HEADER_FILE_INDEX_FILE, this.indexFileName_.get(key)[filePosition]);
    props.put(BioUtilsConstants.INDEX_HEADER_FILE_INDEXED_FILE, this.chromFileName_.get(key)[filePosition]);
    props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_RETENTION_TIME, this.retentionTimeFileName_.get(key)[filePosition]);
    int highestLevel = this.highestMsLevelHash_.get(key).get(filePosition);
    if (msmsInSeveralFiles_) highestLevel = msmsSeveralFilesHighestLevel_;
    props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_MS_LEVEL, String.valueOf(highestLevel));
    if (this.msms_ && highestLevel>1){
      if (msmsInSeveralFiles_){
        props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_MSMS_TYPE, ChromatogramReader.CHROMATOGRAM_HEADER_FILE_MSMS_TYPE_FULL);
      } else
        props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_MSMS_TYPE, ChromatogramReader.CHROMATOGRAM_HEADER_FILE_MSMS_TYPE_PRECURSOR);
      for (int i=2; i<=highestLevel;i++){
        props.put(BioUtilsConstants.INDEX_HEADER_FILE_INDEX_FILE+String.valueOf(i), this.indexFileName_.get(key)[filePosition]+String.valueOf(i));
        props.put(BioUtilsConstants.INDEX_HEADER_FILE_INDEXED_FILE+String.valueOf(i), this.chromFileName_.get(key)[filePosition]+String.valueOf(i));
        props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_RETENTION_TIME+String.valueOf(i), this.retentionTimeFileName_.get(key)[filePosition]+String.valueOf(i));
        if (msmsInSeveralFiles_){
          props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_LOWEST_MZ+String.valueOf(i), String.valueOf(this.lowestMzs_.get(i)));
          props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_HIGHEST_MZ+String.valueOf(i), String.valueOf(this.highestMzs_.get(i)));
          props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_NUMBER_SCANS+String.valueOf(i), String.valueOf(this.getCorrectScanCountMsMsSeveralFiles(i, key, polarity)));   
        }else{
          props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_NUMBER_SCANS+String.valueOf(i), String.valueOf(this.numberOfMs2Scans_.get(key).get(filePosition).get(i)));          
        }
      }
    }
    if (polaritySwitching_){
      String polarityString = null;
      if (polarity==CgDefines.POLARITY_POSITIVE)
        polarityString = GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_POSITIVE;
      else if (polarity==CgDefines.POLARITY_NEGATIVE)
        polarityString = GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_NEGATIVE;
      props.put(GlobalConstants.CHROMATOGRAM_HEADER_FILE_POLARITY_SWITCHED, polarityString);
    }
    FileOutputStream stream = new FileOutputStream(fileName);
    props.store(stream, "Header");
    stream.close();

  }
  
  /**
   * returns the lower/upper m/z threshold for translating the data in integer format
   * depending on the iteration number (memory consumption) and thread number
   * @param iterationNumber the current iteration
   * @param threadNumber the current thread for chrom writing
   * @return the lower/upper m/z threshold for translating the data in integer format
   */
  private RangeInteger getLowerUpperThreshold(int iterationNumber, int threadNumber){
    return this.threadBoundaries_.get(iterationNumber).get(threadNumber);
  }
  
  /**
   * generates the chrom directory and the individual file names
   * @param chromDirName the name and path of the directory that locates the files (without the ".chrom" suffix)
   * @param fileName the name of the file to be translated 
   * @param key the original raw file name (if the raw file contains merged data)
   * @param suffix for MSE files, a suffix is added for the level
   */
  private void generateChromDirectoryAndFiles(String chromDirName, String fileName, String key, String suffix){
    if (polaritySwitching_){
      generateChromFileNames(fileName+RawToChromThread.FILE_SUFFIX_POLARITY_POSITIVE, chromDirName+RawToChromThread.FILE_SUFFIX_POLARITY_POSITIVE+".chrom", key, suffix, 0);
      generateChromFileNames(fileName+RawToChromThread.FILE_SUFFIX_POLARITY_NEGATIVE, chromDirName+RawToChromThread.FILE_SUFFIX_POLARITY_NEGATIVE+".chrom", key, suffix, 1);
    }else{
      generateChromFileNames(fileName, chromDirName+".chrom", key, suffix, 0);
    }
  }
  
  /**
   * generates the chrom directory and the individual file names
   * @param fileName the name of the file to be translated
   * @param chromDirPath the chrom directory that locates the files
   * @param key the original raw file name (if the raw file contains merged data)
   * @param suffix for MSE files, a suffix is added for the level
   * @param the file position (for polarity switching, two files are there)
   */
  private void generateChromFileNames(String fileName, String chromDirPath, String key, String suffix, int position){
    File chromDir = new File(chromDirPath);
    chromDir.mkdir();
    String filePath = chromDir.getAbsolutePath()+File.separator+fileName;
    String[] chromFiles = new String[]{null,null};
    String[] indexFiles = new String[]{null,null};
    String[] headerFiles = new String[]{null,null};
    String[] rtFiles = new String[]{null,null};
    if (chromFileName_.containsKey(key)){
      chromFiles = chromFileName_.get(key);
      indexFiles = indexFileName_.get(key);
      headerFiles = headerFileName_.get(key);
      rtFiles = retentionTimeFileName_.get(key);      
    }
    chromFiles[position] = filePath+".chrom"+suffix;
    indexFiles[position] = filePath+".idx"+suffix;
    headerFiles[position] = filePath+".head"+suffix;
    rtFiles[position] = filePath+".rtt"+suffix;
    
    chromFileName_.put(key,chromFiles);
    indexFileName_.put(key,indexFiles);
    headerFileName_.put(key,headerFiles);
    retentionTimeFileName_.put(key,rtFiles);
  }
  
  private String getCorrespondingThreadDirectory(int polarity, int iteration, int threadNr){
    String dir = m_fileName.substring(0,m_fileName.lastIndexOf("."));
    if (polarity==CgDefines.POLARITY_POSITIVE)
      dir += RawToChromThread.FILE_SUFFIX_POLARITY_POSITIVE;
    else if (polarity==CgDefines.POLARITY_NEGATIVE)
      dir += RawToChromThread.FILE_SUFFIX_POLARITY_NEGATIVE;
    dir += ".chrom/";
    int currentPiece = iteration*numberOfThreads_+threadNr;
    if (currentPiece>0){
      dir +=String.valueOf(currentPiece)+"/";
      File dirFile = new File(dir);
      dirFile.mkdir();
    }
    return dir;
  }

  
  public float getLowerThreshold()
  {
    return this.lowerThreshold_;
  }

  public float getUpperThreshold()
  {
    return this.upperThreshold_;
  }
  
  /**
   * inits the storage hashes for the read scans
   */
  private void initStoreHashes(){
    scanHash_ = new Hashtable<String,CgScan[]>();
    headerHash_ = new Hashtable<String,CgScanHeader>();
    polarity_scanCount_ = new Hashtable<String,Hashtable<Integer,Integer>>();
  }

  /**
   * 
   * @return the interfaces of the RawToChromThreads for the Reader
   */
  private AddScan[] getAdders(){
    AddScan[] adders = new AddScan[translators_.size()];
    for (int i=0;i!=this.translators_.size();i++){
      adders[i] = this.translators_.get(i);
    }
    return adders;
  }
  
  /**
   * @param the position of the file
   * @return the chrom directory
   */
  private String getBaseDir (int filePosition){
    if (polaritySwitching_){
      if (filePosition==0)
        return (this.m_fileName.substring(0,this.m_fileName.lastIndexOf("."))+RawToChromThread.FILE_SUFFIX_POLARITY_POSITIVE+".chrom/");
      else if (filePosition==1)
        return (this.m_fileName.substring(0,this.m_fileName.lastIndexOf("."))+RawToChromThread.FILE_SUFFIX_POLARITY_NEGATIVE+".chrom/");
      else
        return null;
    }else
      return this.m_fileName.substring(0,this.m_fileName.lastIndexOf("."))+".chrom/";
  }

  /**
   * appends one file to another one
   * @param input the file input stream
   * @param output the file to be appended
   * @return the current byte count
   * @throws IOException
   */
  public static long append(InputStream input, OutputStream output)
      throws IOException {
    byte[] buffer = new byte[BUFFER_SIZE];
    long count = 0;
    int n = 0;
    while (-1 != (n = input.read(buffer))) {
      output.write(buffer, 0, n);
      count += n;
    }
    return count;
  }

  /**
   * parses an index file and stores the contents in the provided lines and indices hash tables
   * @param inStream the data input stream of the file
   * @param lines hash containing the line number entries; the key is a counter from 0 to the maximum number of entries
   * @param indices hash containing the byte indices; the key is a counter from 0 to the maximum number of entries
   * @param count the current amount of elements
   * @param previousLength the byte count for previous chrom files (used if file is translated in several slices or threads)
   * @return the current amount of entries
   * @throws IOException thrown if there is something wrong with reading from the data stream
   */
  private int readIndexFile(DataInputStream inStream, Hashtable<Integer,Integer> lines, Hashtable<Integer,Long> indices,
      int count, long previousLength) throws IOException{
    int iteration = 0;
    while (inStream.available()>0){
      int currentLine = inStream.readInt();
      long bytesToSkip = inStream.readLong();
      //fill empty spaces at end of each by thread translated file
      if (iteration==0 && count>0){
        while ((lines.get(count-1)+CgDefines.numberOfEntriesForIndex)<currentLine){
          lines.put(count, lines.get(count-1)+CgDefines.numberOfEntriesForIndex);
          indices.put(count, previousLength);
          count++;
        }
      }
      lines.put(new Integer(count), new Integer(currentLine));
      indices.put(new Integer(count), (new Long(bytesToSkip))+previousLength);
      count++;
      iteration++;
    }
    return count;
  }
  
  /**
   * parses an retention time file and stores the contents in the provided retentionTimes file
   * @param rttStream the data input stream of the file
   * @param retentionTimes the hash containing the retention times; the key is the scan number
   * @throws IOException thrown if there is something wrong with reading from the data stream
   */
  private void readRetentionTimeFile(DataInputStream rttStream, Hashtable<Integer,Float> retentionTimes) throws IOException{
    while (rttStream.available()>0){
      int scanNumber = rttStream.readInt();
      float rtt = rttStream.readFloat();
      retentionTimes.put(scanNumber, rtt);
    }
  }

  /**
   * private class for supervising the threads for chrom translation
   * @author Juergen Hartler
   *
   */
  private class ThreadSupervisor extends TimerTask{
    private int msLevel_;
    
    /**
     * the constructor for the thread supervisor
     * @param msLevel for MSE files, the file is stored in several levels (files)
     * @param currentIteration the number of the current iteration
     */
    public ThreadSupervisor(int msLevel, int currentIteration){
      this.msLevel_ = msLevel;
      startThreads(currentIteration);
    }

    public void run()
    { 
      try{
        if (!oneIterationFinished_)
          handleTimerEvent(msLevel_);
      } catch (Exception ex){
        ex.printStackTrace();
        errorString_ = ex.toString();
        oneIterationFinished_ = true;
      }
    }
    
    private void startThreads(int currentIteration){
      for (int i=0; i!=numberOfThreads_;i++){
        RawToChromThread singleThread = translators_.get(i);
        singleThread.start();   
        quantStatus_.put(i, STATUS_CALCULATING);
      }
    }
  }
  
  /**
   * handles one timer event - checks if the chrom writing threads are finished
   * if all threads are finished, the oneIterationFinished_ is set to true
   * @param msLevel msLevel for MSE files, the file is stored in several levels (files)
   */
  private void handleTimerEvent(int msLevel){
    boolean error = false;
    boolean stopThread = false;
    // find out if some of the threads are finished, read the results, and make the thread available again
    for (int i=0;i!=this.translators_.size();i++){
      boolean isError = false;
      if (this.translators_.get(i)==null) isError = true;
      if (this.translators_.get(i)!=null && this.translators_.get(i).finished()==true && quantStatus_.get(i)!=STATUS_FINISHED){
        RawToChromThread singleThread = this.translators_.get(i);
        if (this.translators_.get(i).getErrorString()!=null) isError=true;
        else{
          if (i==0){
            m_scanCount = singleThread.getm_scanCount();
            polarity_scanCount_ = singleThread.getPolarityScanCount();
            scanCounts_.put(msLevel, m_scanCount);
            polarityScanCounts_.put(msLevel, polarity_scanCount_);
            try{
              singleThread.writeRetentionTimeFile(retentionTimeFileName_,polaritySwitching_);
            } catch (IOException iox){
              iox.printStackTrace();
              errorString_ = iox.getMessage();
              error = true;
              stopThread = true;
            }
            this.numberOfMs2Scans_ = singleThread.getNumberOfMs2Scans();
          }
          singleThread.cleanUp();
          System.gc();
          quantStatus_.put(i, STATUS_FINISHED);
        }
      }
      if (isError){
        error = true;
        stopThread = true;
        if (this.translators_.get(i)!=null) this.errorString_ = this.translators_.get(i).getErrorString();
        else this.errorString_ = "There is no translator for unknown reason - contact the developer";
      }
    }
    
    //check if all of the treads are finished
    boolean allFinished = true;
    if (!error && this.translators_.size()>0){
      for (Integer status : quantStatus_.values()){
        if (status!=STATUS_FINISHED){
          allFinished = false;
          break;
        }
      }
      if (allFinished == true) stopThread = true;
    }
    
    if (stopThread){
      if (!error){
      }
      oneIterationFinished_ = true;
    }
  }
  
  /**
   * returns the adequate scan count for this chrom generation (in case of polarity switched data, only the 
   * corresponding polarity is returned)
   * @param key the file key
   * @param polarity the polarity according to CgDefines
   * @return the adequate scan count for this chrom generation (in case of polarity switched data, only the corresponding polarity is returned)
   */
  private int getCorrectScanCount(String key, int polarity){
    if (polarity==CgDefines.POLARITY_POSITIVE || polarity==CgDefines.POLARITY_NEGATIVE)
      return polarity_scanCount_.get(key).get(polarity);
    else
      return m_scanCount.get(key);    
  }

  /**
   * returns the adequate scan count for this chrom generation for the MSn levels (in case of polarity switched data, only the 
   * corresponding polarity is returned)
   * @param msLevel the MSn level
   * @param key file key
   * @param polarity the polarity according to CgDefines
   * @return the adequate scan count for this chrom generation for the MSn levels (in case of polarity switched data, only the corresponding polarity is returned)
   */
  private int getCorrectScanCountMsMsSeveralFiles(int msLevel, String key, int polarity){
    if (polarity==CgDefines.POLARITY_POSITIVE || polarity==CgDefines.POLARITY_NEGATIVE)
      return polarityScanCounts_.get(msLevel).get(key).get(polarity);
    else
      return scanCounts_.get(msLevel).get(key);
  }
  
  /**
   * 
   * @return true when the mzXML file contains polarity switched data
   */
  public boolean isPolaritySwitched(){
    return this.polaritySwitching_;
  }
  
}
