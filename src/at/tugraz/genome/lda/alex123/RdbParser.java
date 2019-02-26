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
package at.tugraz.genome.lda.alex123;

import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

import at.tugraz.genome.lda.alex123.vos.TargetlistEntry;
import at.tugraz.genome.lda.exception.AlexTargetlistParserException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.vos.QuantVO;

/**
 * Parser for Alex123 target lists
 * @author Juergen Hartler
 *
 */
public class RdbParser
{
  
  /** the file to parse*/
  private String fileName_;
  /** the sequence in which the classes occur*/
  private LinkedHashMap<String,Integer> classSequence_;
  /** the sequence of the found analytes*/
  private Hashtable<String,Vector<String>> analyteSequence_;
  /** the results hash -  in the same format as required for the RdbOutputWriter*/
  private Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> results_;
  
  /**
   * constructor to parse information from the Alex123 RDB file
   * @param fileName the absolute path to the RDB file
   */
  public RdbParser (String fileName) {
    this.fileName_ = fileName;
  }
  
  /**
   * initiates the parsing of the RDB information
   * @throws AlexTargetlistParserException
   */
  public void parse() throws AlexTargetlistParserException{
    results_ = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();
    classSequence_ = new LinkedHashMap<String,Integer>();
    analyteSequence_ = new Hashtable<String,Vector<String>>();
    LinkedHashMap<String,LinkedHashMap<String,String>> analyteSequence = new LinkedHashMap<String,LinkedHashMap<String,String>>();
    LineNumberReader reader = null;
    String line;
    
    try{
      reader = new LineNumberReader(new FileReader(fileName_));
      boolean headerLineFound = false;
      
      int detectorColumn = -1;
      int polarityColumn = -1;
      int speciesColumn = -1;
      int molecularSpeciesColumn = -1;
      int classColumn = -1;
      int idColumn = -1;
      int categoryColumn = -1;
      Hashtable<Integer,Integer> precursorColumns = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Integer> activationColumns = new Hashtable<Integer,Integer>();
      int adductColumn = -1;
      int conflictsColumn = -1;
      int dbNumberColumn = -1;
      int ohNumberColumn = -1;

      String detector;
      String polarity;
      int msLevel;
      String species;
      String molecularSpecies;
      String lClass;
      String id;
      String category;
      String precursor;
      Hashtable<Integer,String> precursors;
      String activation;
      Hashtable<Integer,String> activations;
      String adduct;
      String conflicts;
      int dbNumber;
      int ohNumber;
      
      int classNumber = 1;
      Hashtable<String,Hashtable<String,QuantVO>> resultsOfClass;
      Hashtable<String,QuantVO> resultsOfSpecies;

      while ((line = reader.readLine()) != null) {
        line = line.trim();
//        lineNumber++;
        if (headerLineFound){
          detector = null;
          polarity = null;
          msLevel = 1;
          species = null;
          molecularSpecies = null;
          lClass = null;
          id = null;
          category = null;
          precursor = null;
          precursors = new Hashtable<Integer,String>();
          activation = null;
          activations = new Hashtable<Integer,String>();
          adduct = null;
          conflicts = null;
          dbNumber = -1;
          ohNumber = -1;
          
          String[] columns = line.split("\t");
          if (detectorColumn>=0)
            detector = columns[detectorColumn];
          if (polarityColumn>=0)
            polarity = columns[polarityColumn];
          if (speciesColumn>=0)
            species = columns[speciesColumn];
          if (molecularSpeciesColumn>=0)
            molecularSpecies = columns[molecularSpeciesColumn];
          if (classColumn>=0)
            lClass = columns[classColumn];
          if (idColumn>=0)
            id = columns[idColumn];
          if (categoryColumn>=0)
            category = columns[categoryColumn];
          for (Integer precColumn : precursorColumns.keySet()) {
            precursor = columns[precColumn];
            precursors.put(precursorColumns.get(precColumn), precursor);
            if (precursorColumns.get(precColumn)>msLevel)
              msLevel = precursorColumns.get(precColumn);
          }
          for (Integer actColumn : activationColumns.keySet()) {
            activation = columns[actColumn];
            activations.put(activationColumns.get(actColumn), activation);
            if (activationColumns.get(actColumn)>msLevel)
              msLevel = activationColumns.get(actColumn);
          }
          if (adductColumn>=0)
            adduct = columns[adductColumn];
          if (conflictsColumn>=0)
            conflicts= columns[conflictsColumn];
          if (ohNumberColumn>=0)
            ohNumber = Integer.parseInt(columns[ohNumberColumn]);
          if (dbNumberColumn>=0)
            dbNumber = Integer.parseInt(columns[dbNumberColumn]);
          if (lClass==null || species==null || lClass==null || adduct==null)
            continue;
          if (results_.containsKey(lClass) && results_.get(lClass).containsKey(species) && results_.get(lClass).get(species).containsKey(adduct)) {
            TargetlistEntry entry = (TargetlistEntry)results_.get(lClass).get(species).get(adduct);
            if (msLevel>entry.getMsLevel())
              entry.setMsLevel_(msLevel);
            if (molecularSpecies!=null && (entry.getMolecularSpecies()==null || entry.getMolecularSpecies().length()==0))
              entry.setOriginalMolecularSpecies(molecularSpecies);
            if (precursors.containsKey(2) && (entry.getMs2Precursor()==null || entry.getMs2Precursor().length()==0))
              entry.setMs2Precursor(precursors.get(2));
            if (activations.containsKey(2) && (entry.getMs2Activation()==null || entry.getMs2Activation().length()==0))
              entry.setMs2Activation(activations.get(2));
            if (precursors.containsKey(3) && (entry.getMs3Precursor()==null || entry.getMs3Precursor().length()==0))
              entry.setMs3Precursor(precursors.get(3));
            if (activations.containsKey(3) && (entry.getMs3Precursor()==null || entry.getMs3Precursor().length()==0))
              entry.setMs3Activation(activations.get(3));
          }else {
            //TODO: this is not the appropriate object for this purpose, but it holds all the necessary information
            TargetlistEntry entry = new TargetlistEntry(detector, polarity, msLevel, -1d, null, null, species,
              molecularSpecies, lClass, precursors.containsKey(2) ? precursors.get(2) : null,
              activations.containsKey(2) ? activations.get(2) : null, precursors.containsKey(3) ? precursors.get(3) : null,
              activations.containsKey(3) ? activations.get(3) : null, adduct, null, id, category, conflicts,
              -1, -1, dbNumber, ohNumber, null, null, null, -1, -1, -1, null, null);
          
            resultsOfClass = new Hashtable<String,Hashtable<String,QuantVO>>();
            LinkedHashMap<String,String> analSequ;
            if (classSequence_.containsKey(entry.getAnalyteClass())) {
              resultsOfClass = results_.get(entry.getAnalyteClass());
              analSequ = analyteSequence.get(entry.getAnalyteClass());
            }else {
              classSequence_.put(entry.getAnalyteClass(), classNumber);
              classNumber++;
              analSequ = new LinkedHashMap<String,String>();
              analyteSequence.put(entry.getAnalyteClass(), analSequ);
            }
            resultsOfSpecies = new Hashtable<String,QuantVO>();
            if (resultsOfClass.containsKey(species))
              resultsOfSpecies = resultsOfClass.get(species);
            else
              analSequ.put(species.substring(entry.getAnalyteClass().length()+1), species.substring(entry.getAnalyteClass().length()+1));
            resultsOfSpecies.put(adduct, entry);
            resultsOfClass.put(species, resultsOfSpecies);
            results_.put(entry.getAnalyteClass(), resultsOfClass);
          }
          
        } else if (line.contains(RdbOutputWriter.FILENAME_COLUMN) && line.contains(RdbOutputWriter.MACHINE_COLUMN) &&
            line.contains(RdbOutputWriter.EVIDENCE_LEVEL_COLUMN) && line.contains(RdbOutputWriter.SPECIES_COLUMN) && 
            line.contains(RdbOutputWriter.MOLECULAR_SPECIES_COLUMN) && line.contains(RdbOutputWriter.RT_COLUMN) &&
            line.contains(RdbOutputWriter.PEAK_AREA_COLUMN)){
          String[] columns = line.split("\t");
          for (int columnNr=0; columnNr!=columns.length; columnNr++){
            if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.DETECTOR_COLUMN)){
              detectorColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.POLARITY_COLUMN)){
              polarityColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.SPECIES_COLUMN)){
              speciesColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.MOLECULAR_SPECIES_COLUMN)){
              molecularSpeciesColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.LIPID_CLASS_COLUMN)){
              classColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.LIPID_ID_COLUMN)){
              idColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.LIPID_CATEGORY_COLUMN)){
              categoryColumn = columnNr;
            } else if (columns[columnNr].startsWith(RdbOutputWriter.LEVEL_PREFIX) && 
                (columns[columnNr].endsWith(RdbOutputWriter.ACTIVATION_COLUMN_PART) || columns[columnNr].endsWith(RdbOutputWriter.PRECURSOR_COLUMN_PART))) {
              String headerPart = columns[columnNr].substring(RdbOutputWriter.LEVEL_PREFIX.length());
              int level = Integer.parseInt(headerPart.substring(0,headerPart.indexOf(" ")));
              if (columns[columnNr].endsWith(RdbOutputWriter.PRECURSOR_COLUMN_PART))
                precursorColumns.put(columnNr, level);
              else if (columns[columnNr].endsWith(RdbOutputWriter.ACTIVATION_COLUMN_PART))
                activationColumns.put(columnNr, level);
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.ADDUCT_COLUMN)){
              adductColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.CONFLICTS_COLUMN)){
              conflictsColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.DB_INDEX_COLUMN)){
              dbNumberColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(RdbOutputWriter.OH_INDEX_COLUMN)){
              ohNumberColumn = columnNr;
            }
          }
          headerLineFound = true;
        }
      }
      for (String className : analyteSequence.keySet()) {
        analyteSequence_.put(className, new Vector<String>(analyteSequence.get(className).keySet()));
      }
      
    }catch (IOException | HydroxylationEncodingException ex){
      throw new AlexTargetlistParserException(ex);
    }finally{
      try{if (reader!=null)reader.close();}catch(Exception ex){}
    }
  }

  /**
   * returns the analyte classes in the sequence as they appear in the RDB file
   * @return analyte classes in the sequence as they appear in the RDB file
   */
  public LinkedHashMap<String,Integer> getClassSequence()
  {
    return classSequence_;
  }
  
  
  /**
   * returns the analytes in the sequence as they appear in the RDB file (classwise)
   * @return the analytes in the sequence as they appear in the RDB file (classwise)
   */
  public Hashtable<String,Vector<String>> getAnalyteSequence()
  {
    return analyteSequence_;
  }

  /**
   * the information parsed from the RDB file
   * @return the information parsed from the RDB file
   */
  public Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> getTargetlistInfo()
  {
    return results_;
  }
  
  
  
}
