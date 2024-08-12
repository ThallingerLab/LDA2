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

import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Hashtable;


import java.util.Vector;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.alex123.vos.TargetlistEntry;
import at.tugraz.genome.lda.exception.AlexTargetlistParserException;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;

/**
 * Parser for Alex123 target lists
 * @author Juergen Hartler
 *
 */
public class TargetlistParser
{
  
  /** the file name of the target list*/
  private String fileName_;
  /** true for positive ion mode, false for negative ion mode*/
  private boolean polarity_;
  /** lookup for the notation of isotopes in Alex123 notation to the one in LDA notation*/
  private Hashtable<String,String> isoLookup_;
  /** the results; key is the msLevel; then, a vector of targets*/ 
  private Hashtable<Integer,Vector<TargetlistEntry>> results_;

  private final static String DETECTOR_COLUMN = "Detector";
  private final static String POLARITY_COLUMN = "Polarity";
  private final static String MS_LEVEL_COLUMN = "MS dimension";
  private final static String MZ_COLUMN = "Target m/z";
  private final static String FRAGMENT_COLUMN = "Fragment name";
  private final static String STRUCTURE_COLUMN = "Structure information";
  private final static String SPECIES_COLUMN = "Lipid species";
  private final static String MOLECULAR_SPECIES_COLUMN = "Molecular lipid species";
  private final static String LIPID_CLASS_COLUMN = "Lipid class";
  private final static String MS2_PRECURSOR_COLUMN = "MS2 precursor m/z";
  private final static String MS2_ACTIVATION_COLUMN = "MS2 activation";
  private final static String MS3_PRECURSOR_COLUMN = "MS3 precursor m/z";
  private final static String MS3_ACTIVATION_COLUMN = "MS3 activation";
  private final static String ADDUCT_COLUMN = "Adduct";
  private final static String ID_COLUMN = "Lipid ID";
  private final static String CATEGORY_COLUMN = "Lipid category";
  private final static String CONFLICTS_COLUMN = "Conflicts";
  private final static String CHARGE_COLUMN = "Charge";
  private final static String CARBON_NUMBER_COLUMN = "C index from lipid species";
  private final static String DB_NUMBER_COLUMN = "DB index from lipid species";
  private final static String OH_NUMBER_COLUMN = "OH index from lipid species";
  private final static String SUM_COMPOSITION_COLUMN = "Sum composition from lipid species";
  private final static String FORMULA_COLUMN = "Sum formula from lipid species";
  private final static String FRAGMENT_CARBON_NUMBER_COLUMN = "C index from structure information";
  private final static String FRAGMENT_DB_NUMBER_COLUMN = "DB index from structure information";
  private final static String FRAGMENT_OH_NUMBER_COLUMN = "OH index from structure information";
  private final static String FRAGMENT_SUM_COMPOSITION_COLUMN = "Sum composition from structure information";
  private final static String FRAGMENT_FORMULA_COLUMN = "Sum formula from structure information";
  private final static String FRAGMENT_FORMULA_COLUMN_PARTIAL = "Sum formula from structure inf";
  
  
  /**
   * constructor specifying the Alex123 target list file, and for which polarity the targets shall be extracted
   * @param fileName the file name of the target list
   * @param polarity true for positive ion mode, false for negative ion mode
   */
  public TargetlistParser(String fileName, boolean polarity){
    this.fileName_ = fileName;
    this.polarity_ = polarity;
  }
  
  /**
   * starts the parsing of the Alex123 target list file
   * @throws AlexTargetlistParserException if there is something wrong with the target lists
   */
  public void parse() throws AlexTargetlistParserException{
    String line;
    isoLookup_ = new Hashtable<String,String>();
    results_ = new Hashtable<Integer,Vector<TargetlistEntry>>();
    
    ElementConfigParser elementParser = Settings.getElementParser();
    try (LineNumberReader reader = new LineNumberReader(new FileReader(fileName_));)
    {
      this.isoLookup_ = Settings.getAlexIsoLookup();      
      int lineNumber = 0;
      boolean headerLineFound = false;
      
      int detectorColumn = -1;
      int polarityColumn = -1;
      int msLevelColumn = -1;
      int mzColumn = -1;
      int fragmentColumn = -1;
      int structureColumn = -1;
      int speciesColumn = -1;
      int molecularSpeciesColumn = -1;
      int classColumn = -1;
      int ms2PrecursorColumn = -1;
      int ms2ActivationColumn = -1;
      int ms3PrecursorColumn = -1;
      int ms3ActivationColumn = -1;
      int adductColumn = -1;
      int idColumn = -1;
      int categoryColumn = -1;
      int conflictsColumn = -1;
      int chargeColumn = -1;
      int carbonNumberColumn = -1;
      int dbNumberColumn = -1;
      int ohNumberColumn = -1;
      int sumCompositionColumn = -1;
      int formulaColumn = -1;
      int fragmentCarbonNumberColumn = -1;
      int fragmentDbNumberColumn = -1;
      int fragmentOhNumberColumn = -1;
      int fragmentSumCompositionColumn = -1;
      int fragmentFormulaColumn = -1;
      
      while ((line = reader.readLine()) != null) {
        line = line.trim();
        lineNumber++;
        if (headerLineFound){
          String detector = null;
          String polarity = null;
          int msLevel = -1;
          double mz = -1d;
          String fragment = null;
          String structure = null;
          String species = null;
          String molecularSpecies = null;
          String lipidClass = null;
          String ms2Precursor = null;
          String ms2Activation = null;
          String ms3Precursor = null;
          String ms3Activation = null;
          String adduct = null;
          Hashtable<String,Integer> adductCategorized = null;
          String id = null;
          String category = null;
          String conflicts = null;
          int charge = -1;
          int carbonNumber = -1;
          int dbNumber = -1;
          int ohNumber = -1;
          String sumComposition = null;
          String formula = null;
          Hashtable<String,Integer> formulaCategorized = null;
          int fragmentCarbonNumber = -1;
          int fragmentDbNumber = -1;
          int fragmentOhNumber = -1;
          String fragmentSumComposition = null;
          String fragmentFormula = null;
//          Hashtable<String,Integer> fragmentFormulaCategorized = null;
          
          String[] columns = line.split("\t");
          //reads the single cells and assigns them to the corresponding variables
          for (int columnNr=0; columnNr!=columns.length; columnNr++){
            if (columns[columnNr]==null) continue;
            String entry = columns[columnNr].trim();
            if (entry.length()==0) continue;
            
            if (columnNr==detectorColumn){
              detector = entry;
            } else if (columnNr==polarityColumn){
              polarity = entry;
            } else if (columnNr==msLevelColumn){
              if (entry.equalsIgnoreCase("ms")) msLevel=1;
              else{
                msLevel = parseIntegerEntry(entry.substring(2),MS_LEVEL_COLUMN,lineNumber);
              }
            } else if (columnNr==mzColumn){
              mz = parseDoubleEntry(entry,MZ_COLUMN,lineNumber);
            } else if (columnNr==fragmentColumn){
              fragment = entry;
              if (fragment.startsWith("=")) fragment = fragment.substring(1);
            } else if (columnNr==structureColumn){
              structure = entry;
            } else if (columnNr==speciesColumn){
              species = entry;
            } else if (columnNr==molecularSpeciesColumn){
              molecularSpecies = entry;
            } else if (columnNr==classColumn){
              lipidClass = entry;
            } else if (columnNr==ms2PrecursorColumn){
              parseDoubleEntry(entry,MS2_PRECURSOR_COLUMN,lineNumber);
              ms2Precursor = entry;
            } else if (columnNr==ms2ActivationColumn){
              ms2Activation = entry;
            } else if (columnNr==ms3PrecursorColumn){
              parseDoubleEntry(entry,MS3_PRECURSOR_COLUMN,lineNumber);
              ms3Precursor = entry;
            } else if (columnNr==ms3ActivationColumn){
              ms3Activation = entry;
            } else if (columnNr==adductColumn){
              adduct = entry;
              String adductForCategorization = new String (adduct);
              if (adductForCategorization.endsWith("+") || adductForCategorization.endsWith("-"))
                adductForCategorization = adductForCategorization.substring(0,adductForCategorization.length()-1);
              adductForCategorization = adductForCategorization.replaceAll("\\+ ", " ");
              adductForCategorization = adductForCategorization.replaceAll("- ", " ");
              adductCategorized = StaticUtils.categorizeAdduct(adductForCategorization);
//              adduct = parseChemicalFormula(elementParser, adduct, ADDUCT_COLUMN, lineNumber);
            } else if (columnNr==idColumn){
              id = entry;
            } else if (columnNr==categoryColumn){
              category = entry;
            } else if (columnNr==conflictsColumn){
              conflicts = entry;
            } else if (columnNr==chargeColumn){
              charge = Math.abs(parseIntegerEntry(entry,CHARGE_COLUMN,lineNumber));
            } else if (columnNr==carbonNumberColumn){
              carbonNumber = parseIntegerEntry(entry,CARBON_NUMBER_COLUMN,lineNumber);
            } else if (columnNr==dbNumberColumn){
              dbNumber = parseIntegerEntry(entry,DB_NUMBER_COLUMN,lineNumber);
            } else if (columnNr==ohNumberColumn){
              ohNumber = parseIntegerEntry(entry,OH_NUMBER_COLUMN,lineNumber);
            } else if (columnNr==sumCompositionColumn){
              sumComposition = entry;
            } else if (columnNr==formulaColumn){
              formula = entry;
              formulaCategorized = StaticUtils.categorizeFormula(parseChemicalFormula(elementParser, entry, FORMULA_COLUMN, lineNumber));
            } else if (columnNr==fragmentCarbonNumberColumn){
              fragmentCarbonNumber = parseIntegerEntry(entry,FRAGMENT_CARBON_NUMBER_COLUMN,lineNumber);
            } else if (columnNr==fragmentDbNumberColumn){
              fragmentDbNumber = parseIntegerEntry(entry,FRAGMENT_DB_NUMBER_COLUMN,lineNumber);
            } else if (columnNr==fragmentOhNumberColumn){
              fragmentOhNumber = parseIntegerEntry(entry,FRAGMENT_OH_NUMBER_COLUMN,lineNumber);
            } else if (columnNr==fragmentSumCompositionColumn){
              fragmentSumComposition = entry;
            } else if (columnNr==fragmentFormulaColumn){
              fragmentFormula = parseChemicalFormula(elementParser, entry, FRAGMENT_FORMULA_COLUMN, lineNumber);
//              fragmentFormulaCategorized = StaticUtils.categorizeFormula(fragmentFormula);
            }          
          }
          
          if (msLevel<1)
            throw new AlexTargetlistParserException("An entry must contain information in the \""+MS_LEVEL_COLUMN+"\" column. The entry at line "+lineNumber+" in file "+fileName_+" does not.");
          if (!polarity.equalsIgnoreCase("+") && !polarity.equalsIgnoreCase("-"))
            throw new AlexTargetlistParserException("The \""+POLARITY_COLUMN+"\" column must contain \"+\" or \"-\" as entry. The entry at line "+lineNumber+" in file "+fileName_+" does not.");
          if (polarity_ && !polarity.equalsIgnoreCase("+")) continue;
          else if (!polarity_ && !polarity.equalsIgnoreCase("-")) continue;
          if (mz<0)
            throw new AlexTargetlistParserException("The \""+MZ_COLUMN+"\" column must contain a positive number in double format. The entry at line "+lineNumber+" in file "+fileName_+" is not.");
          if (msLevel>1 && fragment==null) 
            throw new AlexTargetlistParserException("The \""+FRAGMENT_COLUMN+"\" column must not be empty if the \""+MS_LEVEL_COLUMN+"\" is higher than 1. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (msLevel>1 && fragment==null) 
            throw new AlexTargetlistParserException("The \""+STRUCTURE_COLUMN+"\" column must not be empty if the \""+MS_LEVEL_COLUMN+"\" is higher than 1. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (species==null) 
            throw new AlexTargetlistParserException("The \""+SPECIES_COLUMN+"\" column must not be empty. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (msLevel>1 && molecularSpecies==null) 
            throw new AlexTargetlistParserException("The \""+MOLECULAR_SPECIES_COLUMN+"\" column must not be empty if the \""+MS_LEVEL_COLUMN+"\" is higher than 1. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (lipidClass==null) 
            throw new AlexTargetlistParserException("The \""+LIPID_CLASS_COLUMN+"\" column must not be empty. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (msLevel==2 && ms2Precursor==null) 
            throw new AlexTargetlistParserException("The \""+MS2_PRECURSOR_COLUMN+"\" column must not be empty if the \""+MS_LEVEL_COLUMN+"\" equals 2. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (msLevel==2 && ms2Activation==null) 
            throw new AlexTargetlistParserException("The \""+MS2_ACTIVATION_COLUMN+"\" column must not be empty if the \""+MS_LEVEL_COLUMN+"\" equals 2. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (msLevel==3 && ms3Precursor==null) 
            throw new AlexTargetlistParserException("The \""+MS2_PRECURSOR_COLUMN+"\" column must not be empty if the \""+MS_LEVEL_COLUMN+"\" equals 3. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (msLevel==3 && ms3Activation==null) 
            throw new AlexTargetlistParserException("The \""+MS2_ACTIVATION_COLUMN+"\" column must not be empty if the \""+MS_LEVEL_COLUMN+"\" equals 3. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (adduct==null || adductCategorized==null) 
            throw new AlexTargetlistParserException("The \""+ADDUCT_COLUMN+"\" column must not be empty. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (id==null) 
            throw new AlexTargetlistParserException("The \""+ID_COLUMN+"\" column must not be empty. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (category==null) 
            throw new AlexTargetlistParserException("The \""+CATEGORY_COLUMN+"\" column must not be empty. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          if (charge<1)
            throw new AlexTargetlistParserException("An entry in \""+CATEGORY_COLUMN+"\" must be contain a higher integer value than 0. The entry at line "+lineNumber+" in file "+fileName_+" does not.");
          if (formula==null || formulaCategorized==null) 
            throw new AlexTargetlistParserException("The \""+FORMULA_COLUMN+"\" column must not be empty. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");
          //TODO: this is removed meanwhile
//          if (msLevel>1 && (fragmentFormula==null || fragmentFormulaCategorized==null)) 
//            throw new AlexTargetlistParserException("The \""+FRAGMENT_FORMULA_COLUMN+"\" column must not be empty if the \""+MS_LEVEL_COLUMN+"\" is higher than 1. The entry at line "+lineNumber+" in file "+fileName_+" is empty.");

          for (String element : adductCategorized.keySet()){
            int amount = 0;
            if (formulaCategorized.containsKey(element)) amount = formulaCategorized.get(element);
            amount -= adductCategorized.get(element);
            if (amount!=0)
              formulaCategorized.put(element, amount);
            else if (formulaCategorized.containsKey(element))
              formulaCategorized.remove(element);
          }
          TargetlistEntry entry = new TargetlistEntry(detector, polarity, msLevel, mz, fragment, structure, species,
              molecularSpecies, lipidClass, ms2Precursor, ms2Activation, ms3Precursor, ms3Activation, adduct,
              StaticUtils.getFormulaInHillNotation(adductCategorized,true), id, category, conflicts, charge,
              carbonNumber, dbNumber, ohNumber, sumComposition, StaticUtils.getFormulaInHillNotation(formulaCategorized,true),
              formula, fragmentCarbonNumber, fragmentDbNumber, fragmentOhNumber, fragmentSumComposition, fragmentFormula);
//          if (lipidClass.endsWith(" O-") || lipidClass.endsWith(" P-") || lipidClass.equalsIgnoreCase("Cer")) {
//            System.out.println("LDA: "+entry.getAnalyteClass()+" | "+entry.getAnalyteName()+":"+entry.getDbNumber()+" | "+entry.getMolecularSpecies());
//            System.out.println("ALEX: "+entry.getOriginalClassName()+" | "+entry.getSpecies()+" | "+entry.getOriginalMolecularSpecies());
//          }
          Vector<TargetlistEntry> entries = new Vector<TargetlistEntry>();
          if (results_.containsKey(msLevel)) entries = results_.get(msLevel);
          entries.add(entry);
          results_.put(msLevel, entries);
        } else if (line.contains(POLARITY_COLUMN) && line.contains(MS_LEVEL_COLUMN) &&
            line.contains(MZ_COLUMN) && line.contains(LIPID_CLASS_COLUMN) && line.contains(ADDUCT_COLUMN) &&
            line.contains(CHARGE_COLUMN) && line.contains(FORMULA_COLUMN)){
          String[] columns = line.split("\t");
          for (int columnNr=0; columnNr!=columns.length; columnNr++){
            if (columns[columnNr].equalsIgnoreCase(DETECTOR_COLUMN)){
              detectorColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(POLARITY_COLUMN)){
              polarityColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(MS_LEVEL_COLUMN)){
              msLevelColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(MZ_COLUMN)){
              mzColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(FRAGMENT_COLUMN)){
              fragmentColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(STRUCTURE_COLUMN)){
              structureColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(SPECIES_COLUMN)){
              speciesColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(MOLECULAR_SPECIES_COLUMN)){
              molecularSpeciesColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(LIPID_CLASS_COLUMN)){
              classColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(MS2_PRECURSOR_COLUMN)){
              ms2PrecursorColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(MS2_ACTIVATION_COLUMN)){
              ms2ActivationColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(MS3_PRECURSOR_COLUMN)){
              ms3PrecursorColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(MS3_ACTIVATION_COLUMN)){
              ms3ActivationColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(ADDUCT_COLUMN)){
              adductColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(ID_COLUMN)){
              idColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(CATEGORY_COLUMN)){
              categoryColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(CONFLICTS_COLUMN)){
              conflictsColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(CHARGE_COLUMN)){
              chargeColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(CARBON_NUMBER_COLUMN)){
              carbonNumberColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(DB_NUMBER_COLUMN)){
              dbNumberColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(OH_NUMBER_COLUMN)){
              ohNumberColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(SUM_COMPOSITION_COLUMN)){
              sumCompositionColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(FORMULA_COLUMN)){
              formulaColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(FRAGMENT_CARBON_NUMBER_COLUMN)){
              fragmentCarbonNumberColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(FRAGMENT_DB_NUMBER_COLUMN)){
              fragmentDbNumberColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(FRAGMENT_OH_NUMBER_COLUMN)){
              fragmentOhNumberColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(FRAGMENT_SUM_COMPOSITION_COLUMN)){
              fragmentSumCompositionColumn = columnNr;
            } else if (columns[columnNr].equalsIgnoreCase(FRAGMENT_FORMULA_COLUMN) || columns[columnNr].startsWith(FRAGMENT_FORMULA_COLUMN_PARTIAL)){
              fragmentFormulaColumn = columnNr;
            }
          }
          headerLineFound = true;
        }
      }
    }catch (IOException iox){
      throw new AlexTargetlistParserException(iox);
    }catch (ChemicalFormulaException | HydroxylationEncodingException e) {
      throw new AlexTargetlistParserException(e);
    }catch (AlexTargetlistParserException alx){
      throw alx;
    }
  }
  
  /**
   * checks an integer entry for validity and throws an error if there is something wrong
   * @param entry the String entry the putative String entry
   * @param columnName the name of the column for the error message
   * @param lineNumber the line number for the error message
   * @return the entry in integer format
   * @throws AlexTargetlistParserException if the entry is not integer format
   */
  private int parseIntegerEntry(String entry, String columnName, int lineNumber) throws AlexTargetlistParserException{
    try{
      return Integer.parseInt(entry);
    } catch (NumberFormatException nfx){
      throw new AlexTargetlistParserException("Entries in the column \""+columnName+"\" must be integer format. The one at line "+lineNumber+" in file "+fileName_+" is not!");
    }
  }
  
  /**
   * checks an double entry for validity and throws an error if there is something wrong
   * @param entry the String entry the putative String entry
   * @param columnName the name of the column for the error message
   * @param lineNumber the line number for the error message
   * @return the entry in double format
   * @throws AlexTargetlistParserException if the entry is not double format
   */
  private double parseDoubleEntry(String entry, String columnName, int lineNumber) throws AlexTargetlistParserException{
    try{
      return Double.parseDouble(entry);
    } catch (NumberFormatException nfx){
      throw new AlexTargetlistParserException("Entries in the column \""+columnName+"\" must be double format. The one at line "+lineNumber+" in file "+fileName_+" is not!");
    }
  }
  
  /**
   * transfers the Alex123 formula in an LDA readable one (including isotope translation)
   * @param parser object containing the LDA abundances
   * @param formula the Alex123 chemical formula to be parsed
   * @param columnName the name of the column for the error message
   * @param lineNumber the line number for the error message
   * @return the chemical formula in the LDA format
   * @throws AlexTargetlistParserException exception if something is wrong with the entry
   */
  private String parseChemicalFormula(ElementConfigParser parser, String formula, String columnName, int lineNumber) throws AlexTargetlistParserException{
    String formulaString = new String(formula);
    if (formulaString.indexOf("_")!=-1)
      formulaString = formulaString.substring(0,formulaString.indexOf("_"));
    if (formulaString.contains("-")){
      throw new AlexTargetlistParserException("The formula "+formula+" must not contain any negative values! Affected is column \""+columnName+"\" at line "+lineNumber+" in file "+fileName_+".");
    }
    char[] formulaChars = formulaString.toCharArray();
    String formulaToCheck = "";
    boolean isPreviousDigit = false;
    for (int i=0;i!=formulaChars.length;i++){
      char currentChar = formulaChars[i];
      if (isPreviousDigit && !Character.isDigit(currentChar)){
        formulaToCheck+=" ";
      }
      if (currentChar=='['){
        boolean isoFound = false;
        for (String iso : isoLookup_.keySet()){
          if ((i+iso.length()-1)>=formulaChars.length) continue;
          if (formulaString.substring(i,i+iso.length()).equalsIgnoreCase(iso)){
            isoFound = true;
            formulaToCheck += isoLookup_.get(iso);
            i += (iso.length()-1);
            break;
          }
        }
        if (!isoFound)
          throw new AlexTargetlistParserException("The formula \""+formula+"\" contains an istope that is not defined yet. Affected is column \""+columnName+"\" at line "+lineNumber+" in file "+fileName_+".");
      } else
        formulaToCheck+=String.valueOf(currentChar);
      isPreviousDigit = Character.isDigit(currentChar);
      if (!isPreviousDigit && currentChar!=' ' && (i+1)<formulaChars.length && (Character.isUpperCase(formulaChars[i+1]) || formulaChars[i+1]=='[')){
        formulaToCheck += " "; 
      }
    }
    try {
      parser.calculateTheoreticalMass(formulaToCheck, false);
      return formulaToCheck;
    } catch (SpectrummillParserException spx){
      throw new AlexTargetlistParserException(spx.getMessage()+" Affected is column \""+columnName+"\" at line "+lineNumber+" in file "+fileName_+".");
    }
  }


  /**
   * 
   * @return the results; key is the msLevel; then, a vector of targets
   */
  public Hashtable<Integer,Vector<TargetlistEntry>> getResults()
  {
    return results_;
  }
  
  

}
