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
import at.tugraz.genome.vos.MSFinderEntry;
import at.tugraz.genome.vos.MSFinderHitVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MSFinderStructureParser
{
  private String fileName_;
  
  private final static String HEAD_FILE_PATH = "File path";
  private final static String HEAD_FILE_NAME = "File name";
  private final static String HEAD_TITLE = "Title";
  private final static String HEAD_MS1_COUNT = "MS1 count";
  private final static String HEAD_MSMS_COUNT = "MSMS count";
  private final static String HEAD_MZ = "PRECURSORMZ";
  private final static String HEAD_ADDUCT = "PRECURSORTYPE";
  private final static String HEAD_STRUCTURE_RANK = "Structure rank";
  private final static String HEAD_SCORE = "Total score";
  private final static String HEAD_DATABASES = "Databases";
  private final static String HEAD_FORMULA = "Formula";
  private final static String HEAD_ONTOLOGY = "Ontology";
  private final static String HEAD_INCHIKEY = "InChIKey";
  private final static String HEAD_SMILES= "SMILES";
  
  private Vector<MSFinderEntry> results_ = null;

  public MSFinderStructureParser(String fileName){
    this.fileName_ = fileName;
  }
  
  public void parse() throws MSDialException{
    results_ = new Vector<MSFinderEntry>();
    File file = new File(fileName_);
    System.out.println("MS-Finder parsing");
    if (!file.exists()) throw new MSDialException("The file \""+fileName_+"\" does not exist!");
    if (file.isDirectory()) throw new MSDialException("The file \""+fileName_+"\" is a directory!");
    LineNumberReader reader = null;
    try{
      reader = new LineNumberReader(new FileReader(fileName_));
      String line;
      boolean headerFound = false;
      int lineNumber = 0;
      
      int filePathColumn = -1;
      int fileNameColumn = -1;
      int titleColumn = -1;
      int ms1CountColumn = -1;
      int msmsCountColumn = -1;
      int mzColumn = -1;
      int adductColumn = -1;
      //key column; value rank
      Hashtable<Integer,Integer> structureRankColumns = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Integer> scoreColumns = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Integer> databaseColumns = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Integer> formulaColumns = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Integer> ontologyColumns = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Integer> inchikeyColumns = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Integer> smilesColumns = new Hashtable<Integer,Integer>();
      
      String filePath;
      String fileName;
      String title;
      int ms1Count;
      int msmsCount;
      float mz;
      String adduct;
      Hashtable<Integer,String> structureNames;
      Hashtable<Integer,Float> scores;
      Hashtable<Integer,String> databases;
      Hashtable<Integer,String> formulas;
      Hashtable<Integer,String> ontologies;
      Hashtable<Integer,String> inchikeys;
      Hashtable<Integer,String> smiles;

      String structureName;
      float score;
      String database;
      String formula;
      String ontology;
      String inchikey;
      String smile;
      
      while ((line=reader.readLine())!=null){
        lineNumber++;
        int columnCount = 0;
        if (!headerFound){
          if (line.indexOf(HEAD_FILE_PATH)!=-1 && line.indexOf(HEAD_FILE_NAME)!=-1 && line.indexOf(HEAD_TITLE)!=-1 &&
              line.indexOf(HEAD_MS1_COUNT)!=-1 && line.indexOf(HEAD_MSMS_COUNT)!=-1 && line.indexOf(HEAD_MZ)!=-1 &&
              line.indexOf(HEAD_ADDUCT)!=-1 && line.indexOf(HEAD_STRUCTURE_RANK)!=-1 && line.indexOf(HEAD_SCORE)!=-1 &&
              line.indexOf(HEAD_FORMULA)!=-1){
            headerFound = true;
            String[] tokens = line.split("\t");
            int rank = -1;
                        
            for (String header : tokens){
              if (header.equalsIgnoreCase(HEAD_FILE_PATH)) filePathColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_FILE_NAME)) fileNameColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_TITLE)) titleColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_MS1_COUNT)) ms1CountColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_MSMS_COUNT)) msmsCountColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_MZ)) mzColumn = columnCount;
              else if (header.equalsIgnoreCase(HEAD_ADDUCT)) adductColumn = columnCount;
              else if (header.startsWith(HEAD_STRUCTURE_RANK)) {
                rank = Integer.parseInt(header.substring(HEAD_STRUCTURE_RANK.length()).trim());
                structureRankColumns.put(columnCount,rank);
              } else if (header.startsWith(HEAD_SCORE) && rank>0)
                scoreColumns.put(columnCount, rank);
              else if (header.startsWith(HEAD_DATABASES) && rank>0)
                databaseColumns.put(columnCount, rank);
              else if (header.startsWith(HEAD_FORMULA) && rank>0)
                formulaColumns.put(columnCount, rank);
              else if (header.startsWith(HEAD_ONTOLOGY) && rank>0)
                ontologyColumns.put(columnCount, rank);
              else if (header.startsWith(HEAD_INCHIKEY) && rank>0)
                inchikeyColumns.put(columnCount, rank);
              else if (header.startsWith(HEAD_SMILES) && rank>0)
                smilesColumns.put(columnCount, rank);
              columnCount++;
            }
          }
        } else {
          filePath = null;
          fileName = null;
          title = null;
          ms1Count = -1;
          msmsCount = -1;
          mz = -1f;
          adduct = null;
          structureNames = new Hashtable<Integer,String>();
          scores = new Hashtable<Integer,Float>();
          databases = new Hashtable<Integer,String>();
          formulas = new Hashtable<Integer,String>();
          ontologies = new Hashtable<Integer,String>();
          inchikeys = new Hashtable<Integer,String>();
          smiles = new Hashtable<Integer,String>();
          
          String[] tokens = line.split("\t");
          for (String entry : tokens){
            if (filePathColumn>-1 && columnCount==filePathColumn)
              filePath = entry;
            else if (fileNameColumn>-1 && columnCount==fileNameColumn)
              fileName = entry;
            else if (titleColumn>-1 && columnCount==titleColumn)
              title = entry;
            else if (ms1CountColumn>-1 && columnCount==ms1CountColumn)
              ms1Count = Integer.parseInt(entry);
            else if (msmsCountColumn>-1 && columnCount==msmsCountColumn)
              msmsCount = Integer.parseInt(entry);
            else if (mzColumn>-1 && columnCount==mzColumn)
              mz = Float.parseFloat(entry);
            else if (adductColumn>-1 && columnCount==adductColumn)
              adduct = entry;
            else if (structureRankColumns.containsKey(columnCount) && entry!=null && entry.length()>0)
              structureNames.put(structureRankColumns.get(columnCount),entry);
            else if (scoreColumns.containsKey(columnCount) && entry!=null && entry.length()>0)
              scores.put(scoreColumns.get(columnCount),Float.parseFloat(entry));
            else if (databaseColumns.containsKey(columnCount) && entry!=null && entry.length()>0)
              databases.put(databaseColumns.get(columnCount),entry);
            else if (formulaColumns.containsKey(columnCount) && entry!=null && entry.length()>0)
              formulas.put(formulaColumns.get(columnCount),entry);
            else if (ontologyColumns.containsKey(columnCount) && entry!=null && entry.length()>0)
              ontologies.put(ontologyColumns.get(columnCount),entry);
            else if (inchikeyColumns.containsKey(columnCount) && entry!=null && entry.length()>0)
              inchikeys.put(inchikeyColumns.get(columnCount),entry);
            else if (smilesColumns.containsKey(columnCount) && entry!=null && entry.length()>0)
              smiles.put(smilesColumns.get(columnCount),entry);            
            columnCount++;
          }
          if (ms1Count>0 && msmsCount>0 && mz>0 && adduct!=null && structureNames.size()>0) {
            MSFinderEntry vo = new MSFinderEntry(filePath,title, ms1Count, msmsCount, mz, adduct);
            for (int i=0; i!=structureNames.size(); i++) {
              int rank = i+1;
              structureName = null;
              score = -1f;
              database = null;
              formula = null;
              ontology = null;
              inchikey = null;
              smile = null;

              if (structureNames.containsKey(rank)) {
                structureName = structureNames.get(rank);
                if (scores.containsKey(rank) && scores.get(rank)>=0)
                  score = scores.get(rank);
                else throw new MSDialException("There is no score available for the found structure "+structureName+" at line number "+lineNumber);
                if (databases.containsKey(rank) && databases.get(rank)!=null && databases.get(rank).length()>0)
                  database = databases.get(rank);
                if (formulas.containsKey(rank) && formulas.get(rank)!=null && formulas.get(rank).length()>0)
                  formula = formulas.get(rank);
                else throw new MSDialException("There is no formula available for the found structure "+structureName+" at line number "+lineNumber);
                if (ontologies.containsKey(rank) && ontologies.get(rank)!=null && ontologies.get(rank).length()>0)
                  ontology = ontologies.get(rank);
                if (inchikeys.containsKey(rank) && inchikeys.get(rank)!=null && inchikeys.get(rank).length()>0)
                  inchikey = inchikeys.get(rank);
                if (smiles.containsKey(rank) && smiles.get(rank)!=null && smiles.get(rank).length()>0)
                  smile = smiles.get(rank);
                
                vo.addHit(new MSFinderHitVO(rank, structureName, score, database, formula, ontology, inchikey, smile));
              }
              
            }
            if (vo.getHits().size()>0)
              this.results_.add(vo);
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

  
  public Vector<MSFinderEntry> getResults()
  {
    return results_;
  }
  
  
  
}
