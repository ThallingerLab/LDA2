/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.alex123;

import java.io.File;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

import at.tugraz.genome.lda.alex123.vos.TargetlistEntry;
import at.tugraz.genome.lda.exception.AlexTargetlistParserException;

/**
 * Parser for directories containing Alex123 target lists
 * @author Juergen Hartler
 *
 */
public class TargetlistDirParser
{
  /** name of the directory containing the Alex123 target lists*/
  private String dirName_;
  /** was the data acquired in positive ion mode*/
  private boolean positiveIonMode_;
  /** the resulting targets sorted in the same style as in typical LDA target lists
   * first key is the lipid class; the second key is the analyte name; the third key is the modification name */
  private LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>> sortedEntries_;
  
  /**
   * constructor specifying the directory containing the Alex123 target lists and whether data were acquired in positive ion mode or not
   * @param dirName name of the directory containing the Alex123 target lists
   * @param positiveIonMode was the data acquired in positive ion mode
   */
  public TargetlistDirParser(String dirName, boolean positiveIonMode){
    dirName_ = dirName;
    positiveIonMode_ = positiveIonMode;
  }
  
  /**
   * start the parsing of the Alex123 files
   * @throws AlexTargetlistParserException if there is something wrong with the target lists
   */
  public void parse() throws AlexTargetlistParserException {
    File dir = new File(dirName_);
    sortedEntries_ = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>>();
    if (!dir.isDirectory()) throw new AlexTargetlistParserException("The select Alex123 directory is not a directory!");
    File[] files = dir.listFiles();
    Vector<Hashtable<Integer,Vector<TargetlistEntry>>> parsedEntries = new Vector<Hashtable<Integer,Vector<TargetlistEntry>>>();
    for (int i=0; i!=files.length; i++){
      String fileName = files[i].getAbsolutePath();
      if (!fileName.endsWith(".txt")) continue;
      try {
        TargetlistParser parser = new TargetlistParser(fileName,positiveIonMode_);
        parser.parse();
        Hashtable<Integer,Vector<TargetlistEntry>> results = parser.getResults();
        if (results.size()>0) parsedEntries.add(results);
      } catch (AlexTargetlistParserException e) {
        e.printStackTrace();
      }
    }
    //the first key is the lipid class; the second key is the analyte name; the third key is the modification name
    sortedEntries_ = TargetlistDirParser.sortEntriesForLDA(parsedEntries);
  }
  
  /**
   * groups the resulting targets in the same style as in typical LDA target lists
   * @param parsedEntries the entries returned from the TargetlistParser
   * @return the resulting targets sorted in the same style as in typical LDA target lists; first key is the lipid class; the second key is the analyte name; the third key is the modification name
   */
  public static LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>> sortEntriesForLDA(Vector<Hashtable<Integer,Vector<TargetlistEntry>>> parsedEntries){
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>> sortedEntries = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>>();
    Hashtable<String,Boolean> hasClassMsnFragments = new Hashtable<String,Boolean>();
    for (int msLevel=1; msLevel!=5; msLevel++){
      for (Hashtable<Integer,Vector<TargetlistEntry>> oneFile : parsedEntries){
        if (!oneFile.containsKey(msLevel)) continue;
        for (TargetlistEntry entry : oneFile.get(msLevel)){
          if (msLevel==1){
            LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>> classHits = new LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>();
            if (sortedEntries.containsKey(entry.getAnalyteClass())) classHits = sortedEntries.get(entry.getAnalyteClass());
            LinkedHashMap<String,TargetlistEntry> analyteHits = new LinkedHashMap<String,TargetlistEntry>();
            if (classHits.containsKey(entry.getSpecies())) analyteHits = classHits.get(entry.getSpecies());
            if (analyteHits.containsKey(entry.getModName())){
              System.out.println("The species "+entry.getSpecies()+" with the adduct \""+entry.getModName()+"\" is present more than once!");
              continue;
            }
            analyteHits.put(entry.getModName(), entry);
            classHits.put(entry.getSpecies(), analyteHits);
            sortedEntries.put(entry.getAnalyteClass(), classHits);
          //TODO: this is not tested yet with msLevel>2
          }else{
            if (!sortedEntries.containsKey(entry.getAnalyteClass())){
              System.out.println("The Alex123 target list contains the fragment \""+entry.getFragment()+"\" belonging to the species \""+entry.getSpecies()+"\" that belongs to the class \""+entry.getAnalyteClass()+"\" which does not exist in the MS1 entries.");
              continue;
            }
            if (!sortedEntries.get(entry.getAnalyteClass()).containsKey(entry.getSpecies())){
              System.out.println("The Alex123 target list contains the fragment \""+entry.getFragment()+"\" belonging to the species \""+entry.getSpecies()+"\" that  which does not exist in the MS1 entries.");
              continue;
            }
            if (!sortedEntries.get(entry.getAnalyteClass()).get(entry.getSpecies()).containsKey(entry.getModName())){
              System.out.println("The Alex123 target list contains the adduct \""+entry.getFragment()+"\" belonging to the species \""+entry.getSpecies()+"\" belonging to the adduct \""+entry.getModName()+"\" which does not exist in the MS1 entries.");
              continue;
            }
            TargetlistEntry ms1Entry = sortedEntries.get(entry.getAnalyteClass()).get(entry.getSpecies()).get(entry.getModName());
            ms1Entry.addFragment(entry);
            if (!hasClassMsnFragments.containsKey(ms1Entry.getAnalyteClass()))
              hasClassMsnFragments.put(ms1Entry.getAnalyteClass(), true);
          }
        }
      }
    }
    for (String className : sortedEntries.keySet()){
      if (!hasClassMsnFragments.containsKey(className) || !hasClassMsnFragments.get(className)) continue;
      for (LinkedHashMap<String,TargetlistEntry> analyteHits : sortedEntries.get(className).values()){
        for (TargetlistEntry entry : analyteHits.values()){
          entry.setAlex123FragmentsForClass(true);
        }
      }
    }
    return sortedEntries;
  }
  
  /**
   * @return the resulting targets sorted in the same style as in typical LDA target lists; first key is the lipid class; the second key is the analyte name; the third key is the modification name
   */
  public LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>> getResults(){
    return sortedEntries_;
  }
}
