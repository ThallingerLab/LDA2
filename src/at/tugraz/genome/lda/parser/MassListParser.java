/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2024 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.math.BigDecimal;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.alex123.TargetlistDirParser;
import at.tugraz.genome.lda.alex123.TargetlistParser;
import at.tugraz.genome.lda.alex123.vos.TargetlistEntry;
import at.tugraz.genome.lda.exception.AlexTargetlistParserException;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.NoRuleException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.target.export.TargetListExporter;
import at.tugraz.genome.lda.utils.RangeInteger;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class MassListParser
{
	private String quantFile_;
	private float minusTime_;
	private float plusTime_;
	private int amountOfIsotopes_;
	private int isotopesMustMatch_;
	private boolean searchUnknownTime_;
	private float rtShift_;
	private float lowestRetTime_;
	private float highestRetTime_;
	private boolean respectMassShift_;
	private boolean positiveIonMode_;
	
	private LinkedHashMap<String,Integer> classSequence_;
	private LinkedHashMap<String,Vector<String>> analyteSequence_;
	private Hashtable<String,Boolean> adductInsensitiveRtFilter_;
	private Hashtable<String,Boolean> bestMatchBySpectrumCoverage_;
	private Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects_;
	private Hashtable<String,Float> fixedStartTime_;
	private Hashtable<String,Float> fixedStopTime_;
	private Hashtable<String,Integer> ohNumber_;
	private Hashtable<String,RangeInteger> ohRange_;
	
	/**
   * Constructor that parses an LDA mass list file
   * @param quantFile File path of the Quant File
   * @return a vector containing class sequence, analyte sequence, adduct insensitive retention time filter and quantVO objects
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws ExcelInputFileException
   * @throws ChemicalFormulaException
   * @throws RulesException
   * @throws HydroxylationEncodingException
   */
	public MassListParser(String quantFile) throws IOException,SpectrummillParserException,ExcelInputFileException, ChemicalFormulaException, RulesException, HydroxylationEncodingException
	{
		this(quantFile, 0f, 0f, 0, 0, true, 0f, 0f, 0f, true);
	}
	
	/**
   * Constructor that parses an LDA mass list file
   * @param quantFile File path of the Quant File
   * @param minusTime Retention time before tolerance
   * @param plusTime Retention time after tolerance
   * @param amountOfIsotopes Number of isotopes that shall be quantified
   * @param isotopesMustMatch Number of isotopes that must match the theoretical distribution
   * @param searchUnknownTime Search unknown retention time
   * @param rtShift Retention time shift
   * @param lowestRetTime Lowest retention time in the chrom files
   * @param highestRetTime Highest retention time in the chrom files
   * @return a vector containing class sequence, analyte sequence, adduct insensitive retention time filter and quantVO objects
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws ExcelInputFileException
   * @throws ChemicalFormulaException
   * @throws RulesException
   * @throws HydroxylationEncodingException
   */
	public MassListParser(String quantFile, float minusTime, float plusTime,
			int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime,
			float rtShift, float lowestRetTime,
			float highestRetTime) throws IOException,SpectrummillParserException,ExcelInputFileException, ChemicalFormulaException, RulesException, HydroxylationEncodingException
	{
		this(quantFile, minusTime, plusTime, amountOfIsotopes, isotopesMustMatch, searchUnknownTime, rtShift, lowestRetTime, highestRetTime, true);
	}
	
	/**
   * Constructor that parses an LDA mass list file
   * @param quantFile File path of the Quant File
   * @param minusTime Retention time before tolerance
   * @param plusTime Retention time after tolerance
   * @param amountOfIsotopes Number of isotopes that shall be quantified
   * @param isotopesMustMatch Number of isotopes that must match the theoretical distribution
   * @param searchUnknownTime Search unknown retention time
   * @param rtShift Retention time shift
   * @param lowestRetTime Lowest retention time in the chrom files
   * @param highestRetTime Highest retention time in the chrom files
   * @param respectMassShift Take a mass shift range into account
   * @return a vector containing class sequence, analyte sequence, adduct insensitive retention time filter and quantVO objects
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws ExcelInputFileException
   * @throws ChemicalFormulaException
   * @throws RulesException
   * @throws HydroxylationEncodingException
   */
	public MassListParser(String quantFile, float minusTime, float plusTime,
			int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime,
			float rtShift, float lowestRetTime, float highestRetTime, 
			boolean respectMassShift) throws IOException,SpectrummillParserException,ExcelInputFileException, ChemicalFormulaException, RulesException, HydroxylationEncodingException
	{
		this.quantFile_ = quantFile;
		this.minusTime_ = minusTime;
		this.plusTime_ = plusTime;
		this.amountOfIsotopes_ = amountOfIsotopes;
		this.isotopesMustMatch_ = isotopesMustMatch;
		this.searchUnknownTime_ = searchUnknownTime;
		this.rtShift_ = rtShift;
		this.lowestRetTime_ = lowestRetTime;
		this.highestRetTime_ = highestRetTime;
		this.respectMassShift_ = respectMassShift;
		this.classSequence_ = new LinkedHashMap<String,Integer>();
		this.analyteSequence_ = new LinkedHashMap<String,Vector<String>>();
		this.adductInsensitiveRtFilter_ = new Hashtable<String,Boolean>();
		this.bestMatchBySpectrumCoverage_ = new Hashtable<String,Boolean>();
		this.quantObjects_ = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();
		this.fixedStartTime_ = new Hashtable<String,Float>();
		this.fixedStopTime_ = new Hashtable<String,Float>();
		this.ohNumber_ = new Hashtable<String,Integer>();
		this.ohRange_ = new Hashtable<String,RangeInteger>();
		this.parseTargetListFile();
	}
	
	/**
   * Constructor that parses an Alex123 target list file or a directory containing Alex123 target list files
   * @param quantFile the Alex123 file or the directory containing the Alex123 files
   * @param minusTime time tolerance in negative direction
   * @param plusTime time tolerance in positive direction
   * @param amountOfIsotopes the number of isotopes to be quantified
   * @param isotopesMustMatch the number of isotopes that have to matcht the theoretical isotopic distribution 
   * @param searchUnknownTime true if analytes have to be searched where no retention time is present
   * @param basePeakCutoff a relative cutoff value
   * @param rtShift shift of retention times in relation to the entered ones
   * @param lowestRetTime lower hard limit for the retention time
   * @param highestRetTime upper hard limit for the retention time
   * @param respectMassShift take a mass shift range into account
   * @param positiveIonMode should the targets ion positive or in negative ion mode be used for quantitation, only relevant for ALEX123
   * @return a vector containing the targets for quantitation
   * @throws IOException if a file is not there
   * @throws SpectrummillParserException if there is something wrong with the elementconfig.xml
   * @throws AlexTargetlistParserException if there is something wrong with the target lists
   * @throws ChemicalFormulaException if there is something wrong with the chemical formulae
   * @throws RulesException if there is somehting wrong for the rule parsing
   */
	public MassListParser(String quantFile, float minusTime, float plusTime,
			int amountOfIsotopes, int isotopesMustMatch, boolean searchUnknownTime,
			float basePeakCutoff, float rtShift, float lowestRetTime,
			float highestRetTime, boolean respectMassShift, boolean positiveIonMode) throws IOException,SpectrummillParserException, AlexTargetlistParserException, ChemicalFormulaException, RulesException
	{
		this.quantFile_ = quantFile;
		this.minusTime_ = minusTime;
		this.plusTime_ = plusTime;
		this.amountOfIsotopes_ = amountOfIsotopes;
		this.isotopesMustMatch_ = isotopesMustMatch;
		this.searchUnknownTime_ = searchUnknownTime;
		this.rtShift_ = rtShift;
		this.lowestRetTime_ = lowestRetTime;
		this.highestRetTime_ = highestRetTime;
		this.respectMassShift_ = respectMassShift;
		this.positiveIonMode_ = positiveIonMode;
		this.classSequence_ = new LinkedHashMap<String,Integer>();
		this.analyteSequence_ = new LinkedHashMap<String,Vector<String>>();
		this.adductInsensitiveRtFilter_ = new Hashtable<String,Boolean>();
		this.bestMatchBySpectrumCoverage_ = new Hashtable<String,Boolean>();
		this.quantObjects_ = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();
		this.fixedStartTime_ = new Hashtable<String,Float>();
		this.fixedStopTime_ = new Hashtable<String,Float>();
		this.ohNumber_ = new Hashtable<String,Integer>();
		this.ohRange_ = new Hashtable<String,RangeInteger>();
		this.parseAlex123TargetList();
	}

	/**
   * @param quantFile File path of the Quant File
   * @param minusTime Retention time before tolerance
   * @param plusTime Retention time after tolerance
   * @param amountOfIsotopes Number of isotopes that shall be quantified
   * @param isotopesMustMatch Number of isotopes that must match the theoretical distribution
   * @param searchUnknownTime Search unknown retention time
   * @param basePeakCutoff The relative cutoff value in per mille
   * @param rtShift Retention time shift
   * @param lowestRetTime Lowest retention time in the chrom files
   * @param highestRetTime Highest retention time in the chrom files
   * @param respectMassShift Take a mass shift range into account
   * @return a vector containing class sequence, analyte sequence, adduct insensitive retention time filter and quantVO objects
   * @throws IOException
   * @throws SpectrummillParserException
   * @throws ExcelInputFileException
   * @throws ChemicalFormulaException
   * @throws RulesException
   * @throws HydroxylationEncodingException
   */
  private void parseTargetListFile() throws IOException,SpectrummillParserException,ExcelInputFileException, ChemicalFormulaException, RulesException, HydroxylationEncodingException{
    InputStream myxls = new FileInputStream(quantFile_);
    Workbook workbook = null;
    if (quantFile_.endsWith(".xlsx")) workbook = new XSSFWorkbook(myxls);
    else if (quantFile_.endsWith(".xls")) workbook = new HSSFWorkbook(myxls);
    boolean excelOK = false;
    ElementConfigParser aaParser = Settings.getElementParser();
    for (int sheetNumber = 0; sheetNumber!=workbook.getNumberOfSheets(); sheetNumber++){
      Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = new Hashtable<String,Hashtable<String,QuantVO>>();
      Hashtable<String,Hashtable<String,QuantVO>> quantsOfOxClass = new Hashtable<String,Hashtable<String,QuantVO>>();
      
      Vector<String> analytes = new Vector<String>();
      Vector<String> oxAnalytes = new Vector<String>();
      Sheet sheet = workbook.getSheetAt(sheetNumber);
      boolean rtFilterInsensitive = false;
      boolean pickBestMatchBySpectrumCoverage = false;
      int sideChainColumn = -1;
      int doubleBondColumn = -1;
      int molecularSpeciesWithDoubleBondPositionsColumn = -1;
      int oxStateColumn = -1;

      Hashtable<Integer,String> massOfInterestColumns = new Hashtable<Integer,String>();
      Hashtable<String,Hashtable<String,Integer>> adductComposition = new Hashtable<String,Hashtable<String,Integer>>();
      Hashtable<String,Integer> charges = new Hashtable<String,Integer>();
      Hashtable<String,Integer> multi = new Hashtable<String,Integer>();
      int retTimeColumn = -1;
      boolean foundColumns = false;
      float fixedStartTime = 0;
      float fixedEndTime = Float.MAX_VALUE;
      Hashtable<Integer,String> elementColumns = new  Hashtable<Integer,String>();
      int msLevel = 1;
      int ohNumber = LipidomicsConstants.EXCEL_NO_OH_INFO;
      RangeInteger ohRange = null;
      for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
        Row row = sheet.getRow(rowCount);
        String sideChain = "";
        int doubleBonds = -1;
        String oxState = "";
        Hashtable<String,Integer> elementalComposition = new Hashtable<String,Integer>();
        Hashtable <String,Double> massesOfInterest = new Hashtable <String,Double>();
        float retTime = -1;
        String molecularSpeciesWithDB = null;
        Hashtable<Integer,String> possibleElementColumns = new  Hashtable<Integer,String>();
        for (int i=0; row!=null && i!=(row.getLastCellNum()+1);i++){
          Cell cell = row.getCell(i);
          String contents = "";
          Double numeric = null;
          int cellType = -1;
          if (cell!=null) cellType = cell.getCellType();
          if (cellType==Cell.CELL_TYPE_STRING){
            contents = cell.getStringCellValue();
            try{ 
              if (contents!=null)numeric = new Double(contents.replaceAll(",", "."));
            }catch(NumberFormatException nfx){};
          }else if (cellType==Cell.CELL_TYPE_NUMERIC || cellType==Cell.CELL_TYPE_FORMULA){
           numeric = cell.getNumericCellValue();
           contents = String.valueOf(numeric);
          }
          //String contents = sheet.getCell(i,rowCount).getContents();
          if (contents!=null)
            contents = contents.trim();
          if (!foundColumns){
            if (contents.equalsIgnoreCase("Seitenkette")||contents.equalsIgnoreCase("Name")){
              sideChainColumn = i;
            } 
            
            else if (contents.equalsIgnoreCase("dbs")||contents.equalsIgnoreCase("dbs_TAG")){
              doubleBondColumn = i;
            } else if (contents.equalsIgnoreCase(LipidomicsConstants.CHAIN_MOD_COLUMN_NAME)){
              oxStateColumn = i;
            } 

            else if (contents.startsWith("mass")&&contents.contains("(")&&contents.contains(")")){
              String[] formulaAndName = StaticUtils.extractFormulaAndAdductName(contents);
              adductComposition.put(formulaAndName[1],StaticUtils.categorizeFormula(formulaAndName[0]));
              massOfInterestColumns.put(i,formulaAndName[1]);
              charges.put(formulaAndName[1], Integer.parseInt(formulaAndName[2]));
              multi.put(formulaAndName[1], Integer.parseInt(formulaAndName[3]));
            }
            
            else if (contents.equalsIgnoreCase("tR (min)")){
              retTimeColumn = i;
            }
            
            else if (contents.startsWith("Start-RT:")){
              try{
                fixedStartTime = Float.parseFloat(contents.substring("Start-RT:".length()).trim().replaceAll(",", "."));
                fixedStartTime_.put(sheet.getSheetName(), fixedStartTime);
              }catch(NumberFormatException nfx){nfx.printStackTrace();};
            }
            
            else if (contents.startsWith("Stop-RT:")){
              try{
                fixedEndTime = Float.parseFloat(contents.substring("Stop-RT:".length()).trim().replaceAll(",", "."));
                fixedStopTime_.put(sheet.getSheetName(), fixedEndTime);
              }catch(NumberFormatException nfx){};
            }
            
            else if (contents.startsWith("Mass-Trace:")){
              try{
                msLevel = Integer.parseInt(contents.substring("Mass-Trace:".length()).trim().replaceAll(",", "."));
              }catch(NumberFormatException nfx){};  
            }
            
            else if (contents.trim().length()==1||contents.trim().length()==2){
              boolean ok = false;
              if (Character.isUpperCase(contents.trim().toCharArray()[0])){
                if (contents.trim().length()==2){
                  if (Character.isLowerCase(contents.trim().toCharArray()[1]))
                    ok = true;
                }else{
                  ok = true;
                }
                if (ok){
                  String element = contents.trim();
                  if (aaParser.isElementAvailable(element)){
                    possibleElementColumns.put(i, element);
                  }else{
                    if (sideChainColumn>-1)System.out.println("Warning: The elemental column \""+element+"\" does not seem to be a chemical element!");                
                  }
                }
              }
            } else if (contents.trim().equalsIgnoreCase("adductInsensitiveRtFilter")) {
            	rtFilterInsensitive = true;
            } else if (contents.startsWith("OH-Number:")){
              String ohString = contents.substring("OH-Number:".length()).trim().replaceAll(",", ".");
              try{
                ohNumber = Integer.parseInt(ohString);
              //this value need not necessarily to be a number - it might be the hydroxy encoded string
              }catch(NumberFormatException nfx){
                ohNumber = Settings.getLcbHydroxyEncoding().getHydroxyNumber(ohString);
              }
              ohNumber_.put(sheet.getSheetName(), ohNumber);
            } else if (contents.startsWith("OH-Range:")){
              String ohRangeString = contents.substring("OH-Range:".length()).trim().replaceAll(",", ".");
              String[] ohRangeParts = ohRangeString.split("-");
              boolean error = false;
              try{
                if (ohRangeParts.length>2) error = true;
                int start = Integer.parseInt(ohRangeParts[0]);
                int stop = start;
                if (ohRangeParts.length==2)
                  stop = Integer.parseInt(ohRangeParts[1]);
                ohRange = new RangeInteger(start,stop);
              }catch(NumberFormatException nfx){error = true;}
              if (error)
                throw new HydroxylationEncodingException("The value \"OH-Range\" must be a single integer, or a range in the format $lower$-$higher$; the value \""+ohRangeString+"\" in sheet "+sheet.getSheetName()+" does not comply!");
              else
              	ohRange_.put(sheet.getSheetName(), ohRange);
            } else if (contents.trim().equalsIgnoreCase("pickBestMatchBySpectrumCoverage")) {
            	pickBestMatchBySpectrumCoverage = true;
            } else if (contents.equalsIgnoreCase(TargetListExporter.HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS)){
              molecularSpeciesWithDoubleBondPositionsColumn = i;
            }
          }else{            
            if (i==sideChainColumn&&contents!=null&contents.length()>0){
//          for Marlene metabolomics implementation - exclude the if - only "sideChain = contents;" must remain 
              if (numeric!=null){
                sideChain =  String.valueOf((int)Math.round(numeric));
              }else
                sideChain = contents;
            }
            if (i==doubleBondColumn&&contents!=null&&contents.length()>0){
              doubleBonds = numeric.intValue();
            }
            if (i==oxStateColumn&&contents!=null&&contents.length()>0){
              oxState = contents;
            }
            
            if (elementColumns.containsKey(i)&&contents.length()>0){
              int value = 0;
              // this is for columns such as "M" for mass, there the values are float and have
              // to be removed from the chemical formula - however there is no way to figure out
              // directly from Excel if there is an integer entry in the cell or not!
              if (numeric!=null && contents.endsWith(".0")){
                value = (int)Math.round(numeric);
              }else{
                try{
                  value = Integer.parseInt(contents);
                }catch (NumberFormatException nfx3){  
                  System.out.println("Warning: The elemental column \""+elementColumns.get(i)+"\" does not seem to be a chemical element, since float values are there -> I remove it!");
                  elementColumns.remove(i);
                }
              }
              if (value>0)
                elementalComposition.put(elementColumns.get(i), value);
            }

            if (massOfInterestColumns.containsKey(i)&&contents!=null&contents.length()>0){
              double massOfInterest = numeric;
              if (respectMassShift_) massOfInterest += LipidomicsConstants.getMassShift();
              massesOfInterest.put(massOfInterestColumns.get(i), massOfInterest);
            }
            
            if (i==retTimeColumn&&contents!=null&contents.length()>0){
              retTime = numeric.floatValue();
              retTime += rtShift_;
            }
            
            if (i==molecularSpeciesWithDoubleBondPositionsColumn&&contents!=null&contents.length()>0){
              molecularSpeciesWithDB = contents;
            }
          }
        }

        if (sideChainColumn != -1 && possibleElementColumns.size()>0 && 
            massOfInterestColumns.size()>0 && retTimeColumn != -1&&
            !foundColumns){
          foundColumns = true;
          elementColumns = new Hashtable<Integer,String>(possibleElementColumns);
        }
        if (foundColumns&&sideChain!=null&&sideChain.length()>0&&massesOfInterest.size()>0/*&&retTime>0*/){
          float usedMinusTime = new Float(minusTime_);
          if (usedMinusTime<0) usedMinusTime = usedMinusTime*-1f;
          float usedPlusTime = new Float(plusTime_);
          if (usedPlusTime<0) usedPlusTime = usedPlusTime*-1f;
          if (retTime>0||searchUnknownTime_){
            // if a global retention time for a class is set the other retention times have to be recalculated         
            if (fixedStartTime>0 || fixedEndTime<Float.MAX_VALUE && (fixedStartTime<fixedEndTime)){
              if (retTime>0){
                if (fixedStartTime<retTime && retTime<fixedEndTime){
                  if (fixedStartTime>(retTime-usedMinusTime))
                    usedMinusTime = retTime-fixedStartTime;
                  if (fixedEndTime<(retTime+usedPlusTime))
                    usedPlusTime = fixedEndTime-retTime;
                }else{
                // Here I do not know what to do; I guess best is to keep the set RT times, then I can quantify items that are outside
                // the general class range, if necessary!
                }
              }else{
                float startTime = lowestRetTime_;
                float stopTime = highestRetTime_;
                if (fixedStartTime>startTime)
                  startTime = fixedStartTime;
                if (fixedEndTime<stopTime)
                  stopTime = fixedEndTime;
                retTime = (startTime+stopTime)/2;
                usedMinusTime = retTime-startTime;
                usedPlusTime = stopTime-retTime;
              }
            }
            
            int startOh = 0;
            int stopOh = 0;
            if (ohNumber>LipidomicsConstants.EXCEL_NO_OH_INFO) {
              startOh = ohNumber;
              stopOh = ohNumber;
              if (ohRange!=null) {
                startOh = ohRange.getStart();
                stopOh = ohRange.getStop();
              }
            }
            
            //create new quantVO objects for each (ox)modified lipid
            String[] oxStates = oxState.split(";");
            if (oxState.length() == 0) {
              for (int oh=startOh; oh<(stopOh+1); oh++) {
                Hashtable<String,QuantVO> quantsOfAnalyte = new Hashtable<String,QuantVO>();
                String analEncoded = null;
                Hashtable<String,Integer> correctedElementalComposition = new Hashtable<String,Integer>(elementalComposition);
                double ohDiff = 0d;
                int ohToUse = LipidomicsConstants.EXCEL_NO_OH_INFO;
                if (ohNumber>LipidomicsConstants.EXCEL_NO_OH_INFO) {
                  ohToUse = oh;
                  if (oh!=ohNumber) {
                    int oxygens = 0;
                    if (correctedElementalComposition.containsKey("O")) oxygens = correctedElementalComposition.get("O");
                    oxygens += (oh-ohNumber);
                    if (oxygens<0) continue;
                    correctedElementalComposition.put("O", oxygens);
                    ohDiff = (oh-ohNumber)*Settings.getElementParser().getElementDetails("O").getMonoMass();
                  }
                }
                
                for (String modName : massesOfInterest.keySet()){
                  Hashtable<String,Integer> modElements = adductComposition.get(modName);
                  Integer charge = charges.get(modName);
                  double massOfInterest = massesOfInterest.get(modName)+ohDiff/((double)charge);
                  Integer mult = multi.get(modName);
                  
                  String[] formulas = getFormulasAsString(correctedElementalComposition,modElements,mult);
                  String analyteFormula = formulas[0];
                  String modificationFormula = formulas[1];
                  String chemicalFormula = formulas[2];
                  //no negative elements are allowed after an applied modification
                  if (chemicalFormula.indexOf("-")!=-1)
                    continue;
                  Object[] distris = getTheoreticalIsoDistributions(aaParser,isotopesMustMatch_,amountOfIsotopes_,chemicalFormula);
                  @SuppressWarnings("unchecked")
									Vector<Double> mustMatchProbabs = (Vector<Double>)distris[0];
                  @SuppressWarnings("unchecked")
                  Vector<Double> probabs = (Vector<Double>)distris[1];
                  int negativeStartValue = (Integer)distris[2];

                  QuantVO quantVO = new QuantVO(sheet.getSheetName(), sideChain, doubleBonds,
                      ohToUse,analyteFormula, massOfInterest, charge, modName,
                      modificationFormula, retTime, usedMinusTime, usedPlusTime,
                      mustMatchProbabs, probabs,negativeStartValue, "");
                  
                  analEncoded = quantVO.getAnalyteName();
                  quantsOfAnalyte.put(modName, quantVO);
                }
                String analyteName = StaticUtils.generateLipidNameString((analEncoded!=null ? analEncoded : sideChain), doubleBonds, -1, "");
                if (!quantsOfClass.containsKey(analyteName)) {
                  analytes.add(analyteName);
                  quantsOfClass.put(analyteName, quantsOfAnalyte);
                } else if (!hasValidInfoForOmegaAssignment(molecularSpeciesWithDB, retTime)) {
                  System.out.println(String.format("Ignoring duplicate analyte %s in mass list (line %s)!", analyteName, (rowCount+1)));
                }
                if (hasValidInfoForOmegaAssignment(molecularSpeciesWithDB, retTime)) {
                  for (String mod : quantsOfAnalyte.keySet()){
                    try {
                      Vector<FattyAcidVO> chainCombination = StaticUtils.decodeFAsFromHumanReadableName(
                          molecularSpeciesWithDB, Settings.getFaHydroxyEncoding(),Settings.getLcbHydroxyEncoding(), false, null);
                      DoubleBondPositionVO doubleBondPositionVO = new DoubleBondPositionVO(
                          chainCombination, retTime);
                      quantsOfClass.get(analyteName).get(mod).addInfoForOmegaAssignment(doubleBondPositionVO);
                    } catch (LipidCombinameEncodingException ex) {
                      System.out.println(ex.getMessage());
                    }
                  }
                }
              }
            	
            } else {
            for (String oxMod : oxStates)
            {
            	oxMod = oxMod.replaceAll("\\s", "");
            	
	            for (int oh=startOh; oh<(stopOh+1); oh++) {
	              Hashtable<String,QuantVO> quantsOfAnalyte = new Hashtable<String,QuantVO>();
	              Hashtable<String,QuantVO> quantsOfOxAnalyte = new Hashtable<String,QuantVO>();
	              
	              String analEncoded = null;
	              Hashtable<String,Integer> correctedElementalComposition = new Hashtable<String,Integer>(elementalComposition);
	              double ohDiff = 0d;
	              int ohToUse = LipidomicsConstants.EXCEL_NO_OH_INFO;
	              if (ohNumber>LipidomicsConstants.EXCEL_NO_OH_INFO) {
	                ohToUse = oh;
	                if (oh!=ohNumber) {
	                  int oxygens = 0;
	                  if (correctedElementalComposition.containsKey("O")) oxygens = correctedElementalComposition.get("O");
	                  oxygens += (oh-ohNumber);
	                  if (oxygens<0) continue;
	                  correctedElementalComposition.put("O", oxygens);
	                  ohDiff = (oh-ohNumber)*Settings.getElementParser().getElementDetails("O").getMonoMass();
	                }
	              }
	              
	              ModificationParser mp = new ModificationParser(oxMod);
	              mp.parse();
	              correctedElementalComposition = mp.getNewChemicalComposition(correctedElementalComposition);
	              
	              for (String modName : massesOfInterest.keySet()){
	                
	                String analyteClass = sheet.getSheetName();
	                
	                
	                //apply respective (ox)modification to all masses and formulas

                    Hashtable<String,Integer> modElements = adductComposition.get(modName);
                    Integer charge = charges.get(modName);
                    double massOfInterest = mp.getNewMass(massesOfInterest.get(modName)+ohDiff/((double)charge));
                    Integer mult = multi.get(modName);

                    String[] formulas = getFormulasAsString(correctedElementalComposition,modElements,mult);
                    String analyteFormula = formulas[0];
                    String modificationFormula = formulas[1];
                    String chemicalFormula = formulas[2];
                                        
                    //no negative elements are allowed after an applied modification
                    if (chemicalFormula.indexOf("-")!=-1)
                      continue;
                    Object[] distris = getTheoreticalIsoDistributions(aaParser,isotopesMustMatch_,amountOfIsotopes_,chemicalFormula);
                    @SuppressWarnings("unchecked")
                    Vector<Double> mustMatchProbabs = (Vector<Double>)distris[0];
                    @SuppressWarnings("unchecked")
                    Vector<Double> probabs = (Vector<Double>)distris[1];
                    int negativeStartValue = (Integer)distris[2];
                	
                	//add ox-prefix to modified lipids
                	if(!oxMod.equals("")) {
                    	analyteClass = "ox" + analyteClass;
                    }
                	
                	QuantVO quantVO = new QuantVO(analyteClass, sideChain, doubleBonds,
                            ohToUse,analyteFormula, massOfInterest, charge, modName,
                            modificationFormula, retTime, usedMinusTime, usedPlusTime,
                            mustMatchProbabs, probabs,negativeStartValue,oxMod);
                            analEncoded = quantVO.getAnalyteName();
                            
                    /*add oxLipids to another sheet (effectively becoming a tab in the GUI)*/
                            //==
                    if(oxMod.equals("")) {
                    	quantsOfAnalyte.put(modName, quantVO);
                    }else {
                    	quantsOfOxAnalyte.put(modName, quantVO);
                    }
	                	
	              }
	              String analyteName = StaticUtils.generateLipidNameString((analEncoded!=null ? analEncoded : sideChain), doubleBonds, -1, oxMod);
	              
	              /*add oxLipids to another sheet (effectively becoming a tab in the GUI)*/
	              if(oxMod.equals("")) {
	            	  analytes.add(analyteName);
	            	  quantsOfClass.put(analyteName, quantsOfAnalyte);
	              }else {
	            	  oxAnalytes.add(analyteName);
	            	  quantsOfOxClass.put(analyteName, quantsOfOxAnalyte);
	              }              
	            }
            }
          	}
          }
        }
      } 
      quantObjects_.put(sheet.getSheetName(), quantsOfClass);
      analyteSequence_.put(sheet.getSheetName(), analytes);
      classSequence_.put(sheet.getSheetName(),msLevel);
      adductInsensitiveRtFilter_.put(sheet.getSheetName(), rtFilterInsensitive);
      bestMatchBySpectrumCoverage_.put(sheet.getSheetName(), pickBestMatchBySpectrumCoverage);
      
      /*add (if applicable) oxLipids to another sheet (effectively becoming a tab in the GUI)*/
      if(!oxAnalytes.isEmpty())
      {
    	  	quantObjects_.put("ox"+sheet.getSheetName(), quantsOfOxClass);
          analyteSequence_.put("ox"+sheet.getSheetName(), oxAnalytes);
          classSequence_.put("ox"+sheet.getSheetName(),msLevel);
          adductInsensitiveRtFilter_.put("ox"+sheet.getSheetName(), rtFilterInsensitive);
          bestMatchBySpectrumCoverage_.put("ox"+sheet.getSheetName(), pickBestMatchBySpectrumCoverage);
      }
      
      if (foundColumns) excelOK = true;
    }
    myxls.close();
    if (!excelOK) throw new ExcelInputFileException("The Excel file is not valid!");
    checkForIsobaricSpecies(classSequence_,analyteSequence_,quantObjects_);
  }
  
  @Deprecated
  @SuppressWarnings({ "unchecked", "rawtypes" })
  public Vector getResultsVector()
  {
  	Vector results = new Vector();
    results.add(classSequence_);
    results.add(analyteSequence_);
    results.add(adductInsensitiveRtFilter_);
    results.add(bestMatchBySpectrumCoverage_);
    results.add(quantObjects_);
    return results;
  }
  
  /**
   * parses an Alex123 target list file or a directory containing Alex123 target list files
   * @param quantFile the Alex123 file or the directory containing the Alex123 files
   * @param minusTime time tolerance in negative direction
   * @param plusTime time tolerance in positive direction
   * @param amountOfIsotopes the number of isotopes to be quantified
   * @param isotopesMustMatch the number of isotopes that have to matcht the theoretical isotopic distribution 
   * @param searchUnknownTime true if analytes have to be searched where no retention time is present
   * @param basePeakCutoff a relative cutoff value
   * @param rtShift shift of retention times in relation to the entered ones
   * @param lowestRetTime lower hard limit for the retention time
   * @param highestRetTime upper hard limit for the retention time
   * @param positiveIonMode should the targets ion positive or in negative ion mode be used for quantitation
   * @return a vector containing the targets for quantitation
   * @throws IOException if a file is not there
   * @throws SpectrummillParserException if there is something wrong with the elementconfig.xml
   * @throws AlexTargetlistParserException if there is something wrong with the target lists
   * @throws ChemicalFormulaException if there is something wrong with the chemical formulae
   * @throws RulesException if there is somehting wrong for the rule parsing
   */
  private void parseAlex123TargetList() throws IOException,SpectrummillParserException, AlexTargetlistParserException, ChemicalFormulaException, RulesException{
    File file = new File(quantFile_);
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>> sortedEntries = null;
    if (file.isFile() && file.getName().endsWith(".txt")){
      Vector<Hashtable<Integer,Vector<TargetlistEntry>>> parsedEntries = new Vector<Hashtable<Integer,Vector<TargetlistEntry>>>();
      TargetlistParser parser = new TargetlistParser(quantFile_,positiveIonMode_);
      parser.parse();
      parsedEntries.add(parser.getResults());
      sortedEntries = TargetlistDirParser.sortEntriesForLDA(parsedEntries);
    }else if (file.isDirectory()){
      TargetlistDirParser dirParser = new TargetlistDirParser(quantFile_,positiveIonMode_);
      dirParser.parse();
      sortedEntries = dirParser.getResults();
    }
    if (sortedEntries==null || sortedEntries.size()==0)
      throw new AlexTargetlistParserException("There are unusable entries in your target list");

    //now generate the corresponding objects
    ElementConfigParser elementParser = Settings.getElementParser();
    for (String className : sortedEntries.keySet()){
      classSequence_.put(className, 1);
      LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>> classEntries = sortedEntries.get(className);
      Vector<String> analytes = new Vector<String>();
      Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = new Hashtable<String,Hashtable<String,QuantVO>>();
      for (String analyteOriginalName : classEntries.keySet()){
        LinkedHashMap<String,TargetlistEntry> analyteEntries = classEntries.get(analyteOriginalName);
        String sideChain = "";
        int doubleBonds = -1;
        Hashtable<String,QuantVO> quantsOfAnalyte = new Hashtable<String,QuantVO>();
        for (String mod : analyteEntries.keySet()){
          TargetlistEntry entry = analyteEntries.get(mod);
          sideChain = entry.getAnalyteName();
          doubleBonds = entry.getDbs();
          entry.setTimeConstraints(-1, minusTime_, plusTime_);
          
          Hashtable<String,Integer> formAnal = StaticUtils.categorizeFormula(entry.getAnalyteFormula());
          Hashtable<String,Integer> formMod = StaticUtils.categorizeFormula(entry.getModFormula());
          String[] formulas = getFormulasAsString(formAnal,formMod,1);
//          String analyteFormula = formulas[0];
//          String modificationFormula = formulas[1];
          String chemicalFormula = formulas[2];
//          System.out.println(className+StaticUtils.generateLipidNameString(sideChain, doubleBonds)+": "+chemicalFormula);
          Object[] distris = getTheoreticalIsoDistributions(elementParser,isotopesMustMatch_,amountOfIsotopes_,chemicalFormula);
          @SuppressWarnings("unchecked")
					Vector<Double> mustMatchProbabs = (Vector<Double>)distris[0];
          @SuppressWarnings("unchecked")
          Vector<Double> probabs = (Vector<Double>)distris[1];
          int negativeStartValue = (Integer)distris[2];
          entry.setDistributionValues(mustMatchProbabs,probabs,negativeStartValue);
          
          
          quantsOfAnalyte.put(entry.getModName(), entry);
        }
        
        String analyteName = StaticUtils.generateLipidNameString(sideChain, doubleBonds, -1, "");
        analytes.add(analyteName);
        quantsOfClass.put(analyteName, quantsOfAnalyte);
      }
      quantObjects_.put(className, quantsOfClass);
      analyteSequence_.put(className, analytes);
      adductInsensitiveRtFilter_.put(className, false);
      bestMatchBySpectrumCoverage_.put(className, false);
    }
    checkForIsobaricSpecies(classSequence_,analyteSequence_,quantObjects_);
    //TODO: these few lines are only for testing purposes
//    classSequence = new LinkedHashMap<String,Integer>();
//    classSequence.put("TAG", 1);
//    analyteSequence = new Hashtable<String,Vector<String>>();
//    Vector<String> analytesOfClass = new Vector<String>();
//    analytesOfClass.add("40:0");
//    analyteSequence.put("TAG", analytesOfClass);
//    Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects2 = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();
//    Hashtable<String,Hashtable<String,QuantVO>> ofClass = new Hashtable<String,Hashtable<String,QuantVO>>();
//    ofClass.put("40:0", quantObjects.get("TAG").get("40:0"));   
//    quantObjects2.put("TAG", ofClass);
//    results.add(classSequence);
//    results.add(analyteSequence);
//    results.add(adductInsensitiveRtFilter);
//    results.add(quantObjects2);
  }
  
  /**
   * computes chemical formula strings out of the analyte composition and the modification 
   * @param elementalComposition the elemental composition of the analyte
   * @param modElements the elemental composition of the modification
   * @param mult multiplication for dimers, trimers, etc.
   * @return [0] neutral analyte formula; [1] formula of modification; [2] total formula
   */
  private String[] getFormulasAsString(Hashtable<String,Integer> elementalComposition, Hashtable<String,Integer> modElements,
      int mult){
    Hashtable<String,Integer> chemicalFormula = new Hashtable<String,Integer>(elementalComposition);
    for (String element : modElements.keySet()) {
      int amount = modElements.get(element);
      if (chemicalFormula.containsKey(element))
        amount += chemicalFormula.get(element);
      chemicalFormula.put(element, amount);
    }
    String[] formulas = new String[3];
    formulas[0] = StaticUtils.getFormulaInHillNotation(elementalComposition, true);
    formulas[1] = StaticUtils.getFormulaInHillNotation(modElements, true);
    formulas[2] = StaticUtils.getFormulaInHillNotation(chemicalFormula,true);
    return formulas;
  }
  
  /**
   * calculates the positive and the negative isotopic distribution, and decides which on has to be used
   * @param elementParser parser containint the relative abundances of the elements
   * @param isotopesMustMatch number of isotopes that have to fit the theoretical isotopic distribution
   * @param amountOfIsotopes number of isotopes that shall be quantified by the LDA algorithm
   * @param chemicalFormula the chemical formula of the analyte where the isotopic distribution has to match
   * @return [0] Vector<Double> containing the probabilities of the isotopes that must match; [1] Vector<Double> containing probabilities of all isotopes; [2] if the distribution goes in the negative direction - how negative is the lowest isotope 
   * @throws SpectrummillParserException if there is something wrong with the elementconfig.xml
   */
  private Object[] getTheoreticalIsoDistributions(ElementConfigParser elementParser, int isotopesMustMatch, int amountOfIsotopes, String chemicalFormula) throws SpectrummillParserException{  
    boolean negativeDistribution = false;
    Vector<Double> probabs = new Vector<Double>();
    if (amountOfIsotopes<isotopesMustMatch)
      amountOfIsotopes = isotopesMustMatch;
    if (amountOfIsotopes>0){
      Vector<Vector<Double>> bothDistris = elementParser.calculateChemicalFormulaIntensityDistribution(chemicalFormula, amountOfIsotopes+1, false);
      probabs = bothDistris.get(0);
      if (bothDistris.size()>1){
        Vector<Double> negDistri = bothDistris.get(1);
        if (StaticUtils.useNegativeDistribution(probabs,negDistri)){
          probabs = negDistri;
          negativeDistribution = true;
        }
      }

    }else{
      probabs.add(1d);
    }
    Vector<Double> mustMatchProbabs = new Vector<Double>();
    if (isotopesMustMatch>0){
      if (amountOfIsotopes == isotopesMustMatch){
        mustMatchProbabs = new Vector<Double>(probabs);
      }else{
        Vector<Vector<Double>> bothDistris = elementParser.calculateChemicalFormulaIntensityDistribution(chemicalFormula, isotopesMustMatch+1, false);
        mustMatchProbabs = bothDistris.get(0);
        if (bothDistris.size()>1){
          Vector<Double> negDistri = bothDistris.get(1);
          if (StaticUtils.useNegativeDistribution(mustMatchProbabs,negDistri)){
            mustMatchProbabs = negDistri;
            negativeDistribution = true;
          }
        }
      }
    }
    int negativeStartValue = 0;
    if (negativeDistribution){
      negativeStartValue = (mustMatchProbabs.size()*-1)+1;
    }
    
    Object[] distris = new Object[3];
    distris[0] = mustMatchProbabs;
    distris[1] = probabs;
    distris[2] = negativeStartValue;
    return distris;
  }
  
  /**
   * @param molecularSpeciesWithDB String of the potential molecular species with double bond position assignment
   * @param retTime retention time in minutes
   * @return true if a molecular species with double bond position information and a corresponding retention time was found in the mass list
   */
  private boolean hasValidInfoForOmegaAssignment(String molecularSpeciesWithDB, float retTime)
  {
    if (molecularSpeciesWithDB != null && retTime>0) return true;
    return false;
  }
  
  /**
   * checks if the quantitation file contains any isobaric species - if so, this information is stored in the QuantVO
   * @param classSequence the lipid classes in sequence
   * @param analyteSequence the analytes in the lipid classes in sequence
   * @param quantObjects the hash containing the the QuantVO, which pertai information about quantitation
   * @throws RulesException exception thrown if something is wrong with an fragmentation rule
   * @throws IOException if a file is not there
   * @throws SpectrummillParserException if there is something wrong with the elementconfig.xml
   */
  private void checkForIsobaricSpecies(LinkedHashMap<String,Integer> classSequence, LinkedHashMap<String,Vector<String>> analyteSequence,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects) throws RulesException, IOException, SpectrummillParserException{
    Vector<String> classes = new Vector<String>(classSequence.keySet());

    // first, the analytes are put in 10Da bins
    Hashtable<Integer,Vector<QuantVO>> tenDaIsobaricClusters = new Hashtable<Integer,Vector<QuantVO>>();
    for (int i=0; i!=classes.size(); i++){
      String lipidClass = classes.get(i);
      Vector<String> analytes = analyteSequence.get(lipidClass);
      Hashtable<String,Hashtable<String,QuantVO>> quantClass = quantObjects.get(lipidClass);
      for (int j=0; j!=analytes.size(); j++){
        String anal = analytes.get(j);
        Hashtable<String,QuantVO> quantAnal = quantClass.get(anal);
        for (String mod : quantAnal.keySet()){
          QuantVO quant = quantAnal.get(mod);
          float mz = (float)quant.getAnalyteMass();
          int clusterId = getTenDaClusterId(mz);
          Vector<QuantVO> vosOfCluster = new Vector<QuantVO>();
          if (tenDaIsobaricClusters.containsKey(clusterId)) vosOfCluster = tenDaIsobaricClusters.get(clusterId);
          vosOfCluster.add(quant);
          tenDaIsobaricClusters.put(clusterId, vosOfCluster);
        }
      }
    }
    Hashtable<String,String> alreadyUsed = new Hashtable<String,String>(); 
    //second the algorithm looks for analytes which are within the tolerance in the 10Da bins
    for (int i=0; i!=classes.size(); i++){
      String lipidClass = classes.get(i);
      Vector<String> analytes = analyteSequence.get(lipidClass);
      Hashtable<String,Hashtable<String,QuantVO>> quantClass = quantObjects.get(lipidClass);
      for (int j=0; j!=analytes.size(); j++){
        String anal = analytes.get(j);
        Hashtable<String,QuantVO> quantAnal = quantClass.get(anal);
        for (String mod : quantAnal.keySet()){
          QuantVO quant1 = quantAnal.get(mod);
          if (quant1.isQuantifiedByOtherIsobar()) continue;
          
          float tol = LipidomicsConstants.getCoarseChromMzTolerance((float)quant1.getAnalyteMass());
          float sameTol = tol/5f;
          try{
            RulesContainer.getAmountOfChains(StaticUtils.getRuleName(quant1.getAnalyteClass(), quant1.getModName()));
          } catch (NoRuleException nrx) {
            if(j==0) { System.out.println(nrx.getMessage()); }
            continue;
          }
          
          float lowerMz = (float)quant1.getAnalyteMass()-tol;
          float upperMz = (float)quant1.getAnalyteMass()+tol;
          float lowerSame = (float)quant1.getAnalyteMass()-sameTol;
          float upperSame = (float)quant1.getAnalyteMass()+sameTol;
          int clusterId = getTenDaClusterId(lowerMz);
          Vector<QuantVO> toCompare = new Vector<QuantVO>();
          if (tenDaIsobaricClusters.containsKey(clusterId)) toCompare.addAll(tenDaIsobaricClusters.get(clusterId));
          if (clusterId!=getTenDaClusterId(upperMz)){
            clusterId = getTenDaClusterId(upperMz);
            if (tenDaIsobaricClusters.containsKey(clusterId)) toCompare.addAll(tenDaIsobaricClusters.get(clusterId));
          }
          for (QuantVO quant2 : toCompare){
            if (quant1.equals(quant2) || alreadyUsed.containsKey(getUniqueQuantVOString(quant2))) continue;
            
            try{
              RulesContainer.getAmountOfChains(StaticUtils.getRuleName(quant2.getAnalyteClass(), quant2.getModName()));
            } catch (NoRuleException nrx) {
              if(j==0) { System.out.println(nrx.getMessage()); }
              continue;
            }
            
            float mz2 = (float)quant2.getAnalyteMass();
            boolean similarMass = lowerMz<mz2 && mz2<upperMz;
            boolean sameMass = lowerSame<mz2 && mz2<upperSame;
            if (similarMass){
              quant1.addIsobaricSpecies(quant2);
              quant2.addIsobaricSpecies(quant1);
//              if ((quant1.getAnalyteClass().equalsIgnoreCase("P-PE")&&quant1.getIdString().equalsIgnoreCase("36:4"))||
//                  (quant2.getAnalyteClass().equalsIgnoreCase("P-PE")&&quant2.getIdString().equalsIgnoreCase("36:4")))
//                System.out.println(quant1.getAnalyteClass()+quant1.getIdString()+"_"+quant1.getModName()+"/"+quant2.getAnalyteClass()+quant2.getIdString()+"_"+quant2.getModName()+" "+quant1.getAnalyteMass()+"/"+mz2);
              if (sameMass) {
                quant2.setQuantifiedByOtherIsobar(true);
              }
            }
          }
          String idString = getUniqueQuantVOString(quant1);
          alreadyUsed.put(idString, idString);
        }
      }
    }      
  }
  
  /**
   * 
   * @param mz mz value
   * @return a unique integer id for a cluster of the size of 10Da
   */
  private int getTenDaClusterId(float mz){
    return Math.round(Calculator.roundFloat(mz, 0, BigDecimal.ROUND_DOWN));
  }
  
  /**
   * this string is required as unique identifiers, to avoid a repeated comparison
   * @return a unique string separateble from ohter analytes
   */
  private String getUniqueQuantVOString(QuantVO vo){
    return vo.getAnalyteClass()+"_;%"+vo.getIdString()+"_;%"+vo.getModName();
  }

	public LinkedHashMap<String,Integer> getClassSequence()
	{
		return classSequence_;
	}

	public LinkedHashMap<String,Vector<String>> getAnalyteSequence()
	{
		return analyteSequence_;
	}

	public Hashtable<String,Boolean> getAdductInsensitiveRtFilter()
	{
		return adductInsensitiveRtFilter_;
	}

	public Hashtable<String,Boolean> getBestMatchBySpectrumCoverage()
	{
		return bestMatchBySpectrumCoverage_;
	}

	public Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> getQuantObjects()
	{
		return quantObjects_;
	}
	
	public Hashtable<String,Float> getFixedStartTime()
	{
		return fixedStartTime_;
	}
	
	public Hashtable<String,Float> getFixedStopTime()
	{
		return fixedStopTime_;
	}
	
	public Hashtable<String,Integer> getOHNumber()
	{
		return ohNumber_;
	}
	
	public Hashtable<String,RangeInteger> getOHRange()
	{
		return ohRange_;
	}
	
	public boolean isFirstRowRelevant(String cName)
	{
		return getAdductInsensitiveRtFilter().get(cName) != null ||
				getBestMatchBySpectrumCoverage().get(cName) != null ||
				getFixedStartTime().get(cName) != null ||
				getFixedStopTime().get(cName) != null ||
				getOHNumber().get(cName) != null ||
				getOHRange().get(cName) != null;
	}
  
}
