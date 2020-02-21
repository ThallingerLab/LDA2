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

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.lda.exception.MzXMLReaderException;

/**
 * Waters usually stores the data in several mass traces
 * When the data is converted in mzXML, for each mass trace, one mzXML is generated
 * this class merges this mzXML files in one mzXML file
 * @author Juergen Hartler
 *
 */
public class MzXMLMergerForWaters
{
  /** the first mzXML file - other mzXML files end with mzXML2, mzXML3, etc.*/
  private String mzXMLBaseFile_;
  /** the highest stored msLevel; i.e. mzXMLx, e.g., mzXML5*/
  private int msLevels_;
  /** the individual readers for; there is one for each mzXML file*/
  private Hashtable<Integer,MzXMLNextScanReader> readers_;
  /** the filePath to the resulting merged file*/
  private String mergedFileName_;
  
  /**
   * 
   * @param mzXMLBaseFile the name of the first MS level mzXML file
   * @param msLevels  the highest stored msLevel; i.e. mzXMLx, e.g., mzXML5
   */
  public MzXMLMergerForWaters (String mzXMLBaseFile, int msLevels){
    this.mzXMLBaseFile_ = mzXMLBaseFile;
    this.msLevels_ = msLevels;
    mergedFileName_ = mzXMLBaseFile_.substring(0,mzXMLBaseFile_.length()-".mzXML".length())+"_merged.mzXML";
  }
  
  /**
   * conducts the merging
   * @throws IOException when something is wrong with the provided mzXML files (e.g. not there)
   * @throws MzXMLReaderException if something in the mzXML files is not as expected
   */
  public void merge() throws IOException, MzXMLReaderException {
    readers_ = new Hashtable<Integer,MzXMLNextScanReader>();
    // this inits the mzXML readers and lets them read to the first <scan> tag
    int totalNumberOfScans = 0;
    for (int i=1; i<=msLevels_; i++){
      String fileName = new String(mzXMLBaseFile_);
      if (i!=1) fileName += String.valueOf(i);
      File file = new File(fileName);
      if (!file.exists()) continue;
      MzXMLNextScanReader scanReader = new MzXMLNextScanReader(fileName);
      scanReader.openReader(i==1);
      totalNumberOfScans += scanReader.getScanCount();
      readers_.put(i, scanReader);
    }
    BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(mergedFileName_));
    // this writes out the general contents of the mzXML file
    byte[] bytes = replaceScanCount(readers_.get(1).getCachedContents(),totalNumberOfScans).getBytes();
    long currentByteOffset = bytes.length;
    out.write(bytes);
    //this reads the first scans in each mzXML file
    for (int i=1; i<=msLevels_; i++){
      if (!readers_.containsKey(i)) continue;
      readers_.get(i).readNextScan();
    }
    int currentScanNumber = 1;
    String indexEntry = "  <index name=\"scan\">\n";
    String lastScanEntry = null;
    int lastMSLevel = 1;
    Vector<String> closingLines = new Vector<String>();
    boolean foundMS1Scan = false;
    boolean lastIteration = false;
    while (!lastIteration){
      if (allFinished()) lastIteration = true;
      MzXMLNextScanReader scanReader = getReaderAtNextTimePoint();
      if (lastScanEntry!=null){
        indexEntry += createOffsetString(currentScanNumber,currentByteOffset);
        String toWrite = replaceScanNumber(lastScanEntry,currentScanNumber);
        currentScanNumber++;
        if (scanReader==null){
          for (int i=(closingLines.size()-1); i!=-1; i--){
            toWrite += closingLines.get(i);
          }
        } else if (scanReader.getMsLevel()<lastMSLevel){
          for (int i=lastMSLevel; i>scanReader.getMsLevel(); i--){
            toWrite += closingLines.lastElement();
            closingLines.remove(closingLines.size()-1);
          }
        } else if (scanReader.getMsLevel()==lastMSLevel){
          // nothing to do if it is the same MS level
        } else if (scanReader.getMsLevel()>lastMSLevel){
          int indexLastLineBreak = toWrite.substring(0,toWrite.length()-1).lastIndexOf("\n");
          closingLines.add(toWrite.substring(indexLastLineBreak+1));
          toWrite = toWrite.substring(0,indexLastLineBreak+1);
        }
        bytes = toWrite.getBytes();
        currentByteOffset += bytes.length;
        out.write(bytes);
      }
      if (!lastIteration){
        if (foundMS1Scan || scanReader.getMsLevel()==1) {
          foundMS1Scan = true;
          lastScanEntry = scanReader.getCachedContents();
          lastMSLevel = scanReader.getMsLevel();
        }
        scanReader.readNextScan();
      }
    }

    String endTag = "  </msRun>\n";
    bytes = endTag.getBytes();
    currentByteOffset += bytes.length;
    out.write(bytes);
    
    indexEntry += "  </index>\n";
    indexEntry += "  <indexOffset>"+String.valueOf(currentByteOffset)+"</indexOffset>\n";
    out.write(indexEntry.getBytes());
    out.write("</mzXML>\n".getBytes());
    out.close();
    for (int i=1; i<=msLevels_; i++){
      if (!readers_.containsKey(i)) continue;
      readers_.get(i).closeReader();
    }
  }
  
  /**
   * method that stores the bytes to the individual offset
   * @param scanNumber the current scan number
   * @param byteOffset the offset in bytes to this scan
   * @return the offset string
   */
  private String createOffsetString(int scanNumber, long byteOffset){
    return "\t<offset id=\""+String.valueOf(scanNumber)+"\">"+String.valueOf(byteOffset)+"</offset>\n";
  }
  
  /**
   * checks if all of the mzXMLReaders reached the end of the file
   * @return true if all mzXMLReaders reached the end
   */
  private boolean allFinished(){
    boolean allFinished = true;
    for (int i=1; i<=msLevels_; i++){
      if (!readers_.containsKey(i)) continue;
      if (!readers_.get(i).isFinished()){
        allFinished = false;
        break;
      }
    }    
    return allFinished;
  }
  
  /**
   * 
   * @return the mzXMLReader that is the next one in time
   */
  private MzXMLNextScanReader getReaderAtNextTimePoint(){
    MzXMLNextScanReader currentReader = null;
    float lowestTime = Float.MAX_VALUE;
    for (MzXMLNextScanReader reader : readers_.values()){
      if (reader.isFinished()) continue;
      float time = Float.parseFloat(reader.getRetentionTime());
      if (time<lowestTime){
        currentReader = reader;
        lowestTime = time;
      }
    }
    return currentReader;
  }
  
  /**
   * in the merging process, the scans have to be renumbered
   * @param entry one scan entry
   * @param currentScanNumber the new scan number to replace the old one
   * @return the renumbered scan entry
   * @throws MzXMLReaderException an exception if the scan entry contains no "num" attribute
   */
  private String replaceScanNumber(String entry, int currentScanNumber) throws MzXMLReaderException{
    if (entry.indexOf("num=\"")!=-1 || entry.indexOf("num = \"")!=-1){
      String beforeString = "";
      if (entry.indexOf("num=\"")!=-1) beforeString = entry.substring(0,entry.indexOf("num=\"")+"num=\"".length());
      else if (entry.indexOf("num = \"")!=-1) beforeString = entry.substring(0,entry.indexOf("num = \"")+"num = \"".length());
      String restString = entry.substring(beforeString.length());
      restString = restString.substring(restString.indexOf("\""));
      return (beforeString+String.valueOf(currentScanNumber)+restString);
    } else throw new MzXMLReaderException("A mzXML File has no mandatory scan number!");
  }
  
  /**
   * replaces the scanCount with the true number of scans of the merged mzXML
   * @param contents the conventional general contents of the first mzXML
   * @param scanCount the new number of scans
   * @return the corrected conventional general contents of the first mzXML
   */
  private String replaceScanCount(String contents, int scanCount){
    String newContents = "";
    String subContents = new String(contents);
    newContents = subContents.substring(0,subContents.indexOf("scanCount")+"scanCount".length());
    subContents = subContents.substring(subContents.indexOf("scanCount")+"scanCount".length());
    newContents += subContents.substring(0,subContents.indexOf("=")+1);
    subContents = subContents.substring(subContents.indexOf("=")+1);
    newContents += subContents.substring(0,subContents.indexOf("\"")+1);
    subContents = subContents.substring(subContents.indexOf("\"")+1);
    newContents += String.valueOf(scanCount);
    newContents += subContents.substring(subContents.indexOf("\""));
    return newContents;
  }

  /**
   * 
   * @return the filepath to the merged file
   */
  public String getMergedFileName()
  {
    return mergedFileName_;
  }
  
  
}
