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

package at.tugraz.genome.lda.parser;

import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;

import at.tugraz.genome.lda.exception.MzXMLReaderException;

/**
 * This is an mzXML Reader that reads scan entries one by one;
 * The method readNextScan will extract the next scan entry from the mzXML file and
 * provides the original contents and information about the msLevel and the retention time
 * @author Juergen Hartler
 *
 */
public class MzXMLNextScanReader
{
  /** the path to the mzXML file */
  private String fileName_;
  /** the msLevel of the read scan entry */
  private int msLevel_;
  /** the full scan entry contents */
  private String cachedContents_;
  /** the number of scans in the file */
  private int scanCount_;
  /** the retention time of the current scan entry contents*/
  private String retentionTime_;
  /** the last read line */
  private String line_;
  /** was the last scan entry alread read? */
  private boolean finished_;
  /** the line reader for reading lines of the mzXML file */
  private LineNumberReader reader_;
  
  /**
   * constructor setting the file name
   * @param fileName full path to the mzXML file
   */
  public MzXMLNextScanReader(String fileName){
    this.fileName_ = fileName;
  }
  
  /**
   * initiates the line reader and reads until a "<scan" tag is detected
   * @param cacheLinesBeforeScan should the lines before the scan tag be cached (in the first mzXML file this information is used)
   * @throws IOException thrown if something is wrong with the file
   */
  public void openReader(boolean cacheLinesBeforeScan) throws IOException{
    reader_ = new LineNumberReader(new FileReader(fileName_));
    cachedContents_ = "";
    finished_ = false;
    while ((line_ = reader_.readLine()) != null){
      if (line_.trim().startsWith("<scan") || line_.trim().startsWith("< scan")) break;
      //find out the scan count
      if (line_.indexOf("scanCount")!=-1){
        String subLine = line_.trim().substring(line_.trim().indexOf("scanCount")+"scanCount".length());
        subLine = subLine.substring(subLine.indexOf("=")+1);
        subLine = subLine.substring(subLine.indexOf("\"")+1);
        subLine = subLine.substring(0,subLine.indexOf("\"")).trim();
        scanCount_ = Integer.parseInt(subLine);
      }
      if (cacheLinesBeforeScan) cachedContents_ += line_+"\n";
    }
  }
  
  /**
   * reads out the next scan entry
   * @throws IOException thrown if something is wrong with the file
   * @throws MzXMLReaderException thrown if the file scan entry does not contain the required attributes msLevel or retentionTime
   */
  public void readNextScan() throws IOException,MzXMLReaderException{
    msLevel_ = -1;
    cachedContents_ = "";
    retentionTime_ = "";
    if (finished_) return;
    if (line_ !=null && (line_.trim().startsWith("<scan") || line_.trim().startsWith("< scan"))) cachedContents_ += line_+"\n";
    while ((line_ = reader_.readLine()) != null){
      if (line_.indexOf("</msRun>")!=-1 || line_.indexOf("</ msRun>")!=-1 || line_.indexOf("</ msRun >")!=-1){
        finished_ = true;
        return;
      }
      cachedContents_ += line_+"\n";
      if (line_.indexOf("msLevel=\"")!=-1 || line_.indexOf("msLevel = \"")!=-1){
        String msLevelString = "";
        if (line_.indexOf("msLevel=\"")!=-1) msLevelString = line_.substring(line_.indexOf("msLevel=\"")+"msLevel=\"".length());
        else if (line_.indexOf("msLevel = \"")!=-1) msLevelString = line_.substring(line_.indexOf("msLevel = \"")+"msLevel = \"".length());
        msLevelString = msLevelString.substring(0,msLevelString.indexOf("\""));
        try{
          msLevel_ = Integer.parseInt(msLevelString);
        } catch (NumberFormatException nfx){
          throw new MzXMLReaderException("There is something wrong with the msLevel at line "+reader_.getLineNumber()+" with the file: "+fileName_);
        }
      }
      if (line_.indexOf("retentionTime=\"")!=-1 || line_.indexOf("retentionTime = \"")!=-1){
        String retentionTimeString = "";
        if (line_.indexOf("retentionTime=\"")!=-1) retentionTimeString = line_.substring(line_.indexOf("retentionTime=\"")+"retentionTime=\"".length());
        else if (line_.indexOf("retentionTime = \"")!=-1) retentionTimeString = line_.substring(line_.indexOf("retentionTime = \"")+"retentionTime = \"".length());
        retentionTimeString = retentionTimeString.substring(0,retentionTimeString.indexOf("\""));
        char[] rtChars = retentionTimeString.toCharArray();
        for (int i=0; i!=rtChars.length; i++){
          if (Character.isDigit(rtChars[i])){
            retentionTimeString = retentionTimeString.substring(i);
            break;
          }
        }
        rtChars = retentionTimeString.toCharArray();
        for (int i=rtChars.length-1; i!=-1; i--){
          if (Character.isDigit(rtChars[i])){
            retentionTimeString = retentionTimeString.substring(0,i+1);
            break;
          }
        }
        try{
          Double.parseDouble(retentionTimeString);
          retentionTime_ = new String(retentionTimeString);
        } catch (NumberFormatException nfx){
          throw new MzXMLReaderException("There is something wrong with the retention time at line "+reader_.getLineNumber()+" with the file: "+fileName_);
        }
      }
      if (line_.trim().startsWith("</scan") || line_.trim().startsWith("</ scan")) break;
    }
  }
  
  /**
   * 
   * @return the scan entry in its original format
   */
  public String getCachedContents()
  {
    return cachedContents_;
  }
  
  
  /**
   * 
   * @return the amount of scans
   */
  public int getScanCount()
  {
    return scanCount_;
  }

  /**
   * closes the reader and releases the file handle
   */
  public void closeReader(){
    try {if (reader_!=null) reader_.close();}catch (IOException e) {}
  }

  /**
   * 
   * @return the MS level of the current scan entry
   */
  public int getMsLevel()
  {
    return msLevel_;
  }

  /**
   * 
   * @return the retention time of the current scan entry
   */
  public String getRetentionTime()
  {
    return retentionTime_;
  }

  /**
   * 
   * @return true if there are no more scan entries
   */
  public boolean isFinished()
  {
    return finished_;
  }
  
}
