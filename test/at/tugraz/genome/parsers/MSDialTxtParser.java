/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2019 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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
package at.tugraz.genome.parsers;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Hashtable;
import java.util.Vector;

import at.tugraz.genome.exception.MSDialException;
import at.tugraz.genome.vos.MSDialEntry;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MSDialTxtParser
{
  private String fileName_;

  private final static String HEAD_ID = "PeakID";
  private final static String HEAD_TITLE = "Title";
  private final static String HEAD_SCANS = "Scans";
  private final static String HEAD_RT_START = "RT left(min)";
  private final static String HEAD_RT = "RT (min)";
  private final static String HEAD_RT_STOP = "RT right (min)";
  private final static String HEAD_MZ = "Precursor m/z";
  private final static String HEAD_AREA = "Area";
  private final static String HEAD_ADDUCT = "Adduct";
  private final static String HEAD_ISOTOPE = "Isotope";
  private final static String HEAD_SCORE = "Total score";
  private final static String HEAD_SN = "S/N";
  private final static String HEAD_MSMS = "MSMS spectrum";
  
  
  private Vector<MSDialEntry> results_ = null;
  private Vector<MSDialEntry> ms1Only_ = null;
  
  public MSDialTxtParser(String fileName){
    this.fileName_ = fileName;
  }

  public void parse() throws MSDialException{
    results_ = new Vector<MSDialEntry>();
    ms1Only_ = new Vector<MSDialEntry>();
    File file = new File(fileName_);
    System.out.println("MS-Dial parsing");
    if (!file.exists()) throw new MSDialException("The file \""+fileName_+"\" does not exist!");
    if (file.isDirectory()) throw new MSDialException("The file \""+fileName_+"\" is a directory!");
    LineNumberReader reader = null;
    try{
      reader = new LineNumberReader(new FileReader(fileName_));
      String line;
      MSDialEntry vo;
      boolean headerFound = false;
      
      int idColumn = -1;      
      int titleColumn = -1;
      int scansColumn = -1;
      int rtStartColumn = -1;
      int rtColumn = -1;
      int rtStopColumn = -1;
      int mzColumn = -1;
      int areaColumn = -1;
      int adductColumn = -1;
      int isotopeColumn = -1;
      int scoreColumn = -1;
      int signalNoiseColumn = -1;
      int msmsSpectrumColumn  = -1;
      
      String id;
      String name;
      int scans;
      float rtStart;
      float rt;
      float rtStop;
      double mz;
      float area;
      String adduct;
      String isotope;
      float score;
      float signalNoise;
      String msmsSpectrum;
      
      int lineNumber = 0;
      
      
      while ((line=reader.readLine())!=null){
        lineNumber++;
        int columnCount = 0;
        if (!headerFound){
          if (line.indexOf(HEAD_ID)!=-1 && line.indexOf(HEAD_TITLE)!=-1 && line.indexOf(HEAD_RT)!=-1 &&
              line.indexOf(HEAD_MZ)!=-1 && line.indexOf(HEAD_ADDUCT)!=-1 && line.indexOf(HEAD_ISOTOPE)!=-1 &&
              line.indexOf(HEAD_SCORE)!=-1){
            headerFound = true;
            String[] tokens = line.split("\t");
                        
            for (String header : tokens){
              if (header.equalsIgnoreCase(HEAD_ID)) idColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_TITLE)) titleColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_SCANS)) scansColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_RT_START)) rtStartColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_RT)) rtColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_RT_STOP)) rtStopColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_MZ)) mzColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_AREA)) areaColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_ADDUCT)) adductColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_ISOTOPE)) isotopeColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_SCORE)) scoreColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_SN)) signalNoiseColumn = columnCount; 
              else if (header.equalsIgnoreCase(HEAD_MSMS)) msmsSpectrumColumn = columnCount;
              columnCount++;
            }
          }
        } else {
          id = null;
          name = null;
          scans = -1;
          rtStart = -1f;
          rt = -1f;
          rtStop = -1f;
          mz = -1d;
          area = -1f;
          adduct = null;
          isotope = null;
          score = -1f;
          signalNoise = -1f;
          msmsSpectrum = null;
          
          String[] tokens = line.split("\t");
          for (String entry : tokens){
            if (idColumn>-1 && columnCount==idColumn) {
              id = entry;
            }else if (titleColumn>-1 && columnCount==titleColumn) {
              name = entry;
            }else if (scansColumn>-1 && columnCount==scansColumn) {
              scans = Integer.parseInt(entry);
            }else if (rtStartColumn>-1 && columnCount==rtStartColumn) {
              rtStart = Float.parseFloat(entry);
            }else if (rtColumn>-1 && columnCount==rtColumn) {
              rt = Float.parseFloat(entry);
            }else if (rtStopColumn>-1 && columnCount==rtStopColumn) {
              rtStop = Float.parseFloat(entry);
            }else if (mzColumn>-1 && columnCount==mzColumn) {
              mz = Double.parseDouble(entry);
            }else if (areaColumn>-1 && columnCount==areaColumn) {
              area = Float.parseFloat(entry);
            }else if (adductColumn>-1 && columnCount==adductColumn) {
              adduct = entry;
            }else if (isotopeColumn>-1 && columnCount==isotopeColumn) {
              isotope = entry;
            }else if (scoreColumn>-1 && columnCount==scoreColumn) {
              try{
                score = Float.parseFloat(entry);
              } catch(NumberFormatException nfx) {
                System.out.println("ScoreColumn: "+scoreColumn+"; "+lineNumber);
                throw nfx;
              }
            }else if (signalNoiseColumn>-1 && columnCount==signalNoiseColumn) {
              signalNoise = Float.parseFloat(entry);
            }else if (msmsSpectrumColumn>-1 && columnCount==msmsSpectrumColumn) {
              msmsSpectrum = entry;
            }
            columnCount++;
          }
          if (id!= null && name!=null && !name.equalsIgnoreCase("Unknown") /*&& msmsSpectrum!=null && msmsSpectrum.length()>0*/) {
            if (rt<0 || mz<0 || area<0 || adduct==null || isotope==null || score<0) {
              throw new MSDialException ("There is something wrong with the entry on line number "+lineNumber+": "+line);
            }
            //the entry in the adduct column is not reliable;
            adduct = name.substring(name.lastIndexOf(";")+1).trim();
            name = name.substring(0,name.lastIndexOf(";"));
            vo = new MSDialEntry(id, name, scans, rtStart, rt, rtStop, mz, area, adduct, isotope, score, signalNoise, msmsSpectrum);
            if (!vo.getIsotope().equalsIgnoreCase("M + 0")) {
              throw new MSDialException("!!! Identification with "+HEAD_ID+" "+id+" is based on another isotope: "+vo.getIsotope()+" - this hit is discarded !!!");
            }
            if (msmsSpectrum!=null && msmsSpectrum.length()>0)
              results_.add(vo);
            else
              ms1Only_.add(vo);
          }

        }
      }
    } catch (IOException iox){
      throw new MSDialException(iox);
    } finally {
      if (reader!=null){
        try {
          reader.close();
        }
        catch (IOException e){}
      }  
    }
  }

  public Vector<MSDialEntry> getResults()
  {
    return results_;
  }

  public Vector<MSDialEntry> getMs1Only()
  {
    return ms1Only_;
  }
  
  public Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> getStructuredResults(){
    return getStructuredResults(results_);
  }
  
  public Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> getStructuredMS1Only(){
    return getStructuredResults(ms1Only_);
  }  

  
  private static Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> getStructuredResults(Vector<MSDialEntry> results){
    Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> structured = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>>();
    String classString;
    String analyte;
    String mod;
    for (MSDialEntry entry : results) {
      classString = entry.getLdaClassName()!=null ? entry.getLdaClassName() : entry.getDialClassName();
      analyte = entry.getLdaMs1Name()!=null ? entry.getLdaMs1Name() : entry.getDialMs1Name();
      mod = entry.getAdduct();
      Hashtable<String,Hashtable<String,Vector<MSDialEntry>>> structClass = new Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>();
      if (structured.containsKey(classString))
        structClass = structured.get(classString);
      Hashtable<String,Vector<MSDialEntry>> structAnal = new Hashtable<String,Vector<MSDialEntry>>();
      if (structClass.containsKey(analyte))
        structAnal = structClass.get(analyte);
      Vector<MSDialEntry> hits = new Vector<MSDialEntry>();
      if (structAnal.containsKey(mod))
        hits = structAnal.get(mod);
      hits.add(entry);
      structAnal.put(mod, hits);
      structClass.put(analyte, structAnal);
      structured.put(classString, structClass);
    }
    return structured;
  }
  
}
