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

package at.tugraz.genome.lda;

import java.io.File;
import java.util.Timer;
import java.util.TimerTask;

import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JProgressBar;


import at.tugraz.genome.lda.swing.BatchQuantificationTable;
import at.tugraz.genome.lda.swing.BatchQuantificationTableModel;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.RawQuantificationPairVO;
import at.tugraz.genome.maspectras.utils.StringUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class BatchQuantThread extends Thread
{

  private QuantificationThread quantThread_;
  private RawToMzxmlThread rawmzThread_;
  private MzxmlToChromThread mzThread_;
  
  private BatchQuantificationTable quantTable_;
  private BatchQuantificationTableModel quantTableModel_;
  private JProgressBar progressBar_;
  private JLabel quantifyingLabel_;
//  private float mzTolerance_;
  private float minusTime_;
  private float plusTime_;
  private Timer timer_;
  private boolean finished_;
  private int amountOfIsotopes_;
  int isotopesMustMatch_;
  private boolean searchUnknownTime_;
  private float basePeakCutoff_;
  private float rtShift_;
  
  private int currentLine_;
  private boolean readFromRaw_;
  int oneThird_;
  private int numberOfProcessors_;
  private boolean areWiffPresent_;
  private Vector<RawQuantificationPairVO> generatedMzXMLsFromWiff_;
  
  public BatchQuantThread(BatchQuantificationTable quantTable, BatchQuantificationTableModel quantTableModel, 
      JProgressBar progressBar, JLabel quantifyingLabel,//float mzTolerance,
      float minusTime, float plusTime,int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime,
      float basePeakCutoff, float rtShift, int numberOfProcessors) {
    this.quantTable_ = quantTable;
    this.quantTableModel_ = quantTableModel;
    this.progressBar_ = progressBar;
    this.quantifyingLabel_ = quantifyingLabel;
//    this.mzTolerance_ = mzTolerance;
    this.minusTime_ = minusTime;
    this.plusTime_ = plusTime;
    currentLine_ = 0;
    this.quantThread_ = null;
    this.rawmzThread_ = null;
    this.mzThread_ = null;
    this.readFromRaw_ = false;
    this.oneThird_ = 100/(3*quantTableModel.getRowCount());
    this.amountOfIsotopes_ = amountOfIsotopes;
    this.isotopesMustMatch_ = isotopesMustMatch;
    this.searchUnknownTime_ = searchUnknownTime;
    this.basePeakCutoff_ = basePeakCutoff;
    rtShift_ = rtShift;
    numberOfProcessors_ = numberOfProcessors;
    this.generatedMzXMLsFromWiff_ = new Vector<RawQuantificationPairVO>();
    this.areWiffPresent_ = false;
  }
    
  public void run(){
    timer_ = new java.util.Timer();
    timer_.schedule(new ThreadSupervisor(), 10, 1000);
  }
  
  public boolean finished(){
    if (finished_)
      this.timer_.cancel();
    return this.finished_;
  }
  
  private class ThreadSupervisor extends TimerTask{

    public void run()
    {
      handleTimerEvent();
    }   
  }
  
  private void handleTimerEvent(){
    if (this.currentLine_<this.quantTableModel_.getRowCount()){
      if (this.rawmzThread_!=null && this.rawmzThread_.finished()){
        RawQuantificationPairVO filePair = quantTableModel_.getDataByRow((this.currentLine_));
        if (this.rawmzThread_.getErrorString()!=null&&this.rawmzThread_.getErrorString().length()>0){
          filePair.setStatus("ERROR");
          quantTableModel_.fileQuantificationError(filePair);
          this.currentLine_++;
        }else{
          this.rawmzThread_ = null;
          if (filePair.getRawFileName().endsWith(".wiff")){
            areWiffPresent_ = true;
            Vector<File> filesToTranslate = getMzXMLFilesOfWiffConversion(filePair.getRawFile().getAbsolutePath());
            for (File fileToTranslate : filesToTranslate){
              boolean isThere = false;
              for (RawQuantificationPairVO attachedFile : generatedMzXMLsFromWiff_){
                if (fileToTranslate.getAbsolutePath().equalsIgnoreCase(attachedFile.getRawFile().getAbsolutePath()) &&
                    filePair.getQuantFile().getAbsolutePath().equalsIgnoreCase(attachedFile.getQuantFile().getAbsolutePath())) isThere = true;
              }
              if (!isThere) generatedMzXMLsFromWiff_.add(new RawQuantificationPairVO(fileToTranslate,filePair.getQuantFile(),true));
            }
            quantTableModel_.fileQuantified(filePair);
            this.currentLine_ = this.currentLine_+1;
          } else {
            filePair.setStatus("Trans to chrom");
            this.quantTableModel_.fireTableDataChanged();
            this.quantifyingLabel_.setText("Translating "+filePair.getRawFileName()+" to chrom");
            this.progressBar_.setValue((this.currentLine_*100)/this.quantTableModel_.getRowCount()+oneThird_);
            this.readFromRaw_ = true;
            String mzXMLFilePath = filePair.getRawFile().getAbsolutePath().substring(0,filePair.getRawFile().getAbsolutePath().lastIndexOf("."))+".mzXML";
            mzThread_ = new MzxmlToChromThread(mzXMLFilePath);
            mzThread_.start();
          }
        }
      }
      if (this.mzThread_!=null && this.mzThread_.finished()){
        RawQuantificationPairVO filePair = quantTableModel_.getDataByRow((this.currentLine_));
        if (this.mzThread_.getErrorString()!=null&&this.mzThread_.getErrorString().length()>0){
          filePair.setStatus("ERROR");
          quantTableModel_.fileQuantificationError(filePair);
          this.currentLine_++;
        }else{
          filePair.setStatus("Quantifying");
          this.quantifyingLabel_.setText("Quantifying "+filePair.getRawFileName()+" with "+filePair.getQuantFileName());
          this.progressBar_.setValue((this.currentLine_*100)/this.quantTableModel_.getRowCount()+2*oneThird_);
          quantThread_ = new QuantificationThread(filePair.getRawFile().getAbsolutePath(), filePair.getQuantFile().getAbsolutePath(),
              LipidDataAnalyzer.getResultFilePath(filePair.getRawFile().getAbsolutePath(), filePair.getQuantFile().getAbsolutePath()),
              //this.mzTolerance_,
              minusTime_,plusTime_,this.amountOfIsotopes_,this.isotopesMustMatch_,this.searchUnknownTime_,this.basePeakCutoff_,rtShift_, 
              numberOfProcessors_);
          quantThread_.start();
        }
        this.mzThread_ = null;
        if (readFromRaw_ || filePair.isFromWiff()){
          this.readFromRaw_ = false;
          RawToMzxmlThread.deleteMzXMLFiles(filePair.getRawFile().getAbsolutePath().substring(0,filePair.getRawFile().getAbsolutePath().lastIndexOf("."))+".mzXML");
        }
      }
      if (this.quantThread_!=null && !this.quantThread_.finished() && 
          (this.quantThread_.getErrorString()==null||this.quantThread_.getErrorString().length()==0)){   
        if (quantThread_.getTotalAmountOfLipids()>0&&quantThread_.getCurrentLipidCount()>0){
          RawQuantificationPairVO filePair = quantTableModel_.getDataByRow((this.currentLine_));
          this.progressBar_.setValue((this.currentLine_*100)/this.quantTableModel_.getRowCount()+2*oneThird_+((oneThird_*(quantThread_.getCurrentLipidCount()-1))/quantThread_.getTotalAmountOfLipids()));
          filePair.setStatus("Quantifying "+quantThread_.getCurrentLipid()+" ("+quantThread_.getCurrentLipidCount()+"/"+quantThread_.getTotalAmountOfLipids()+")");
          this.quantTable_.repaint();
        }  
      }
      
      if (this.quantThread_!=null && this.quantThread_.finished()){
        RawQuantificationPairVO filePair = quantTableModel_.getDataByRow((this.currentLine_));
        if (this.quantThread_.getErrorString()!=null&&this.quantThread_.getErrorString().length()>0){
          String errorString = quantThread_.getErrorString();
          if (errorString.indexOf(":")>-1) errorString = errorString.substring(errorString.indexOf(":")+1);
          filePair.setStatus("ERROR: "+errorString);
          quantTableModel_.fileQuantificationError(filePair);
        }else{
          filePair.setStatus("Finished");
          quantTableModel_.fileQuantified(filePair);
        }
        this.quantThread_ = null;
        this.quantTable_.getSelectionModel().removeSelectionInterval(0, this.quantTable_.getRowCount());
        this.progressBar_.setValue((this.currentLine_*100)/this.quantTableModel_.getRowCount());
        this.currentLine_ = this.currentLine_+1;
      }
      if (this.currentLine_<this.quantTableModel_.getRowCount()&& this.rawmzThread_ == null && this.mzThread_==null && this.quantThread_ == null){
        readFromRaw_ = false;
        boolean threadStarted = false;
        RawQuantificationPairVO filePair = quantTableModel_.getDataByRow((this.currentLine_));
        String suffix = filePair.getRawFile().getAbsolutePath().substring(filePair.getRawFile().getAbsolutePath().lastIndexOf(".")+1);
        if (suffix.equalsIgnoreCase("RAW")||suffix.equalsIgnoreCase("d")||suffix.equalsIgnoreCase("wiff")){
          File rawFile = new File(filePair.getRawFile().getAbsolutePath());
          if ((rawFile.isFile()&& ((Settings.getReadWPath()!=null&&Settings.getReadWPath().length()>0 && !LipidomicsConstants.isMS2()) ||
               (Settings.getMsConvertPath()!=null && Settings.getMsConvertPath().length()>0))) ||
               (rawFile.isDirectory()&& ((Settings.getMassWolfPath()!=null&&Settings.getMassWolfPath().length()>0 && !LipidomicsConstants.isMS2()) ||
                   (Settings.getMassPlusPlusPath()!=null && Settings.getMassPlusPlusPath().length()>0)))){
            File headerFile = new File(StringUtils.getChromFilePaths(filePair.getRawFile().getAbsolutePath())[1]);
            File mzXMLFile = new File(filePair.getRawFile().getAbsolutePath().substring(0,filePair.getRawFile().getAbsolutePath().length()-suffix.length())+"mzXML");
            if (!headerFile.exists()&&!mzXMLFile.exists()){
              boolean isMassPlusPlus = false;
              filePair.setStatus("Trans to mzXML");
              this.quantTable_.scrollToCenter(this.currentLine_, 0);
              this.quantTable_.getSelectionModel().setSelectionInterval(this.currentLine_, this.currentLine_);
              this.quantTable_.repaint();
              this.quantifyingLabel_.setText("Translating "+filePair.getRawFileName()+" to mzXML");
              this.progressBar_.setValue((this.currentLine_*100)/this.quantTableModel_.getRowCount());
//              String execCommand = "";
//              if (rawFile.isFile())
//                execCommand =  "\""+LipidomicsConstants.getReadWPath()+"\" \""+ filePair.getRawFile().getAbsolutePath()+"\" p";
//              if (rawFile.isDirectory())
//                execCommand =  "\""+LipidomicsConstants.getMassWolfPath()+"\" --mzXML \""+ filePair.getRawFile().getAbsolutePath()+"\" \""+mzXMLFile+"\"";
              String[] params = new String[3];
              if (rawFile.isFile()){
                if (Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0){
                  params = BatchQuantThread.getMsConvertParams(filePair.getRawFile().getAbsolutePath());
                } else if (Settings.getReadWPath()!=null&&Settings.getReadWPath().length()>0){
                  params[0] = Settings.getReadWPath();
                  params[1] = filePair.getRawFile().getAbsolutePath();
                  params[2] = "p";
                }  
              }
              if (rawFile.isDirectory()){
                if (suffix.equalsIgnoreCase("RAW")){
                  if (Settings.getMassPlusPlusPath()!=null&&Settings.getMassPlusPlusPath().length()>0){
                    params = new String[8];
                    params[0] = Settings.getMassPlusPlusPath();
                    params[1] = "-in";
                    params[2] = filePair.getRawFile().getAbsolutePath();
                    params[3] = "-out";
                    params[4] = "mzxml";
                    params[5] = mzXMLFile.getAbsolutePath();
                    params[6] = "-sample";
                    params[7] = "0";
                    if (LipidomicsConstants.isMS2()) isMassPlusPlus = true;
                  }else if (Settings.getMassWolfPath()!=null&&Settings.getMassWolfPath().length()>0){
                    params = new String[4];
                    params[0] = Settings.getMassWolfPath();
                    params[1] = "--mzXML";
                    params[2] = filePair.getRawFile().getAbsolutePath();
                    params[3] = mzXMLFile.getAbsolutePath();  
                  }
                }else if(suffix.equalsIgnoreCase("d")){
                  if (Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0){
                    params =BatchQuantThread.getMsConvertParams(filePair.getRawFile().getAbsolutePath());
                  }                    
                }
              }  
              this.rawmzThread_ = new RawToMzxmlThread(params,isMassPlusPlus);
              rawmzThread_.start();
              threadStarted = true;
            }
          }else{
            filePair.setStatus("ERROR");
            quantTableModel_.fileQuantificationError(filePair);
            this.currentLine_++;
            threadStarted=true;
            return;            
          }
        }
//        if (LipidomicsConstants.getReadWPath()==null||LipidomicsConstants.getReadWPath().length()<1){
//          filePair.setStatus("ERROR");
//          quantTableModel_.fileQuantificationError(filePair);
//          this.currentLine_++;
//          threadStarted=true;
//          return;
//        }
        if (suffix.equalsIgnoreCase("mzXML")){
          File headerFile = new File(StringUtils.getChromFilePaths(filePair.getRawFile().getAbsolutePath())[1]);
          if (!headerFile.exists()){
            filePair.setStatus("Trans to chrom");
            this.quantTable_.scrollToCenter(this.currentLine_, 0);
            this.quantTable_.getSelectionModel().setSelectionInterval(this.currentLine_, this.currentLine_);
            this.quantTable_.repaint();
            this.quantifyingLabel_.setText("Translating "+filePair.getRawFileName()+" to chrom");
            this.progressBar_.setValue((this.currentLine_*100)/this.quantTableModel_.getRowCount()+oneThird_);
            mzThread_ = new MzxmlToChromThread(filePair.getRawFile().getAbsolutePath());
            mzThread_.start();
            threadStarted = true;
          }
        }
        if (!threadStarted){
          filePair.setStatus("Quantifying");
          this.quantTable_.scrollToCenter(this.currentLine_, 0);
          this.quantTable_.getSelectionModel().setSelectionInterval(this.currentLine_, this.currentLine_);
          this.quantTable_.repaint();
          this.quantifyingLabel_.setText("Quantifying "+filePair.getRawFileName()+" with "+filePair.getQuantFileName());
          this.progressBar_.setValue((this.currentLine_*100)/this.quantTableModel_.getRowCount()+2*oneThird_);
          quantThread_ = new QuantificationThread(filePair.getRawFile().getAbsolutePath(), filePair.getQuantFile().getAbsolutePath(),
              LipidDataAnalyzer.getResultFilePath(filePair.getRawFile().getAbsolutePath(), filePair.getQuantFile().getAbsolutePath()),
              //this.mzTolerance_,
              minusTime_,plusTime_,this.amountOfIsotopes_,this.isotopesMustMatch_,this.searchUnknownTime_,this.basePeakCutoff_,rtShift_,
              numberOfProcessors_);
          quantThread_.start();
          threadStarted = true;
        }
        if (!threadStarted){
          filePair.setStatus("ERROR");
          quantTableModel_.fileQuantificationError(filePair);
          this.currentLine_++;
          
        }
      }
    }
    if (this.currentLine_>=this.quantTableModel_.getRowCount()){
      if (areWiffPresent_){
        areWiffPresent_ = false;
        quantTableModel_.clearFiles();
        quantTableModel_.addFiles(this.generatedMzXMLsFromWiff_);
        this.oneThird_ = 100/(3*quantTableModel_.getRowCount());
        currentLine_ = 0;
      }else this.finished_ = true;
    }
  }
  
  
  public static String[] getMsConvertParams(String fileToTranslate){
    String[] params = new String[5];
    params[0] = Settings.getMsConvertPath();
    params[1] = "--mzXML";
    params[2] = fileToTranslate;
    params[3] = "-o";
    params[4] = StaticUtils.extractDirName(fileToTranslate);                  
    return params;
  }
  
  /**
   * 
   * @param filePath the path to the wiff file
   * @return list of mzXML files that were generated from this wiff file
   */
  public static Vector<File> getMzXMLFilesOfWiffConversion(String filePath){
    Vector<File> filesToTranslate = new Vector<File>();
    String fileStart = filePath.substring(0,filePath.lastIndexOf("."));
    File[] filesOfDir = new File(StringUtils.extractDirName(fileStart)).listFiles();
    for (int i=0; i!=filesOfDir.length; i++){          
      if (filesOfDir[i].getAbsolutePath().startsWith(fileStart)&&filesOfDir[i].getAbsolutePath().endsWith(".mzXML")){
        filesToTranslate.add(filesOfDir[i]);
      }
    }
    return filesToTranslate;
  }

}