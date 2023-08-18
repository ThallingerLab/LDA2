/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2020 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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


package at.tugraz.genome.lda.export;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFFont;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFFont;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.analysis.ComparativeAnalysis;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.export.vos.AnalyteOmegaInfoVO;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.DoubleIntegerVO;
import at.tugraz.genome.lda.vos.IsotopicLabelVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * Class for exporting the omega mass lists
 * @author Juergen Hartler
 *
 */
public class OmegaMasslistExporter extends LDAExporter
{
  
  /** the position of the columns*/
  private final static int COLUMN_NAME = 0;
  private final static int COLUMN_COLON = 1;
  private final static int COLUMN_DBS = 2;
  private final static int COLUMN_MOL = 3;
  private final static int COLUMN_FIRST_ELEMENT = 4;
  /** the file name to export the mass lists*/
  private String fileName_;
  
  private final static String HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS = "mol. species";
  
  /**
   * constructor setting the mass list file
   * @param fileName the file name of the exported omega mass list
   */
  public OmegaMasslistExporter(String fileName) {
    this.fileName_ = fileName;
  }
  
  public static String getMolecularSpeciesWithDoubleBondPositionsHeader() {
    return HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS;
  }
  
  /**
   * exports the omega mass list to an Excel file
   * @param classSequence the sequence of the lipid classes
   * @param correctAnalyteSequence the sequence of the analytes within a class; key: lipid class
   * @param quantObjects must not be null: the original quantification objects - these one are required for all the other species that have not been labeled
   * @param acceptedMolecules the actual molecules that have been accepted in the heat map, and are actually labeled species; first key: analyte class; second key: unlabeled analyte; third key: id of the label (an analyte can be assigned to several labels); the labels that were found for this analyte and this labelID (can be different RT and double labelings)
   * @param analysisModule the module containing all information extracted from the result Excel files
   * @param sameOmegaLabels hash grouping labels that have the same omega position; first key: omega position; second key: label identifier; value: label information
   * @throws ExportException when something goes wrong with the export
   * @throws LipidCombinameEncodingException when there is something wrong with the lipid name encoding
   */
  public void export(LinkedHashMap<String,Integer> classSequence, Hashtable<String,Vector<String>> correctAnalyteSequence, Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects,
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<AnalyteOmegaInfoVO>>>> acceptedMolecules, ComparativeAnalysis analysisModule, 
      Hashtable<Integer,Hashtable<String,IsotopicLabelVO>> sameOmegaLabels) throws ExportException, LipidCombinameEncodingException{
    //TODO: this is for checking whether the correct molecules enter this area
//    for (String lClass : acceptedMolecules.keySet()) {
//      for (String molecule : acceptedMolecules.get(lClass).keySet()) {
//        for (String label : acceptedMolecules.get(lClass).get(molecule).keySet()) {
//          for (AnalyteOmegaInfoVO info : acceptedMolecules.get(lClass).get(molecule).get(label)) {
//            //if (info.isSingelDoubleBond())
//            //System.out.println(lClass+" "+molecule+": "+info.getLabelDetection()+"     "+label+" ; "+info.getAppliedLabelId());
//          }
//        }
//        //System.out.println("Accepted: "+lClass+" "+molecule);
//      }
//    }
    double stdev = getMeanStandardDeviationOfSingleDoubleBonds(acceptedMolecules, analysisModule);
    System.out.println("stdev: "+stdev);
    BufferedOutputStream out = null;
    XSSFWorkbook workbook = null;
    try{
      ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
      parser.parse();
      
      out = new BufferedOutputStream(new FileOutputStream(fileName_));
      workbook = new XSSFWorkbook();
      XSSFCellStyle headerStyle = getHeaderStyle(workbook);   
      for (String cName : classSequence.keySet()) {
        Sheet sheet = workbook.createSheet(cName);
        writeMassListForSheet(sheet, headerStyle, correctAnalyteSequence.get(cName), quantObjects.get(cName),acceptedMolecules.get(cName), analysisModule, cName,
           stdev, analysisModule.getNrOfChainsOfClass().get(cName), sameOmegaLabels);
        
      }
      
    } catch (SpectrummillParserException | FileNotFoundException | ChemicalFormulaException e) {
      throw new ExportException(e);
    } finally {
      if (workbook!=null) {
        try {
          workbook.write(out);
          workbook.close();
        } catch (IOException e) {throw new ExportException(e);}
      }
      if (out!=null)
        try {out.close();} catch (IOException e) {throw new ExportException(e);}
    }
  }
  
  
  /**
   * writes an Excel sheet containing the mass list of a single analyte class
   * @param sheet Excel sheet to write in
   * @param headerStyle the style of the header column
   * @param analytes the sequence of the analytes
   * @param quantObjects the values read from a 'conventional' mass list
   * @param labeled the detected species carrying a label; first key: unlabeled analyte; second key: id of the label (an analyte can be assigned to several labels); the labels that were found for this analyte and this labelID (can be different RT and double labelings)
   * @param analysisModule the module containing all information extracted from the result Excel files
   * @param className analyte class name
   * @param stdev mean standard deviation of the species with a single double bond
   * @param nrChains the number of chains this analyte class does have
   * @param sameOmegaLabels hash grouping labels that have the same omega position; first key: omega position; second key: label identifier; value: label information
   * @throws ChemicalFormulaException thrown when there is something wrong with the chemical formula
   * @throws ExportException when something goes wrong with the export
   * @throws LipidCombinameEncodingException when there is something wrong with the lipid name encoding
   */
  @SuppressWarnings("unchecked")
  private void writeMassListForSheet(Sheet sheet, XSSFCellStyle headerStyle, Vector<String> analytes, Hashtable<String,Hashtable<String,QuantVO>> quantObjects,
      Hashtable<String,Hashtable<String,Vector<AnalyteOmegaInfoVO>>> labeled, ComparativeAnalysis analysisModule, String className, 
      double stdev, int nrChains, Hashtable<Integer,Hashtable<String,IsotopicLabelVO>> sameOmegaLabels) throws ChemicalFormulaException, ExportException, LipidCombinameEncodingException {
    //first key: encoded omega chain information; second key: RT
    Hashtable<String,Hashtable<Double,DoubleIntegerVO>> sumSpeciesHash;
    //first key: chain idents; second key: expName; third key: RT
    Hashtable<String,Hashtable<String,Hashtable<Double,DoubleIntegerVO>>> molSpeciesHash;
    ResultAreaVO areaVO;
    Hashtable<String,Integer> elementsLabeled;
    int nr;
    String chainsWithoutLabel;
    String labeledChainsOnly;
    String encodedOmegaChainInformation;
    Double rt;
    DoubleIntegerVO countRts;
    double rtTolerance = stdev*3d;
    double labelTol = stdev*2d;
    Hashtable<String,Hashtable<String,List<Double>>> rtsOfOneSpecies;
    Hashtable<String,Hashtable<String,List<Double>>> rtsSpeciesAndOmega;
    
    Vector<Object> elModCharge = getAvailableElementsAndModificationsPlusCharge(quantObjects);
    Vector<String> elements = (Vector<String>)elModCharge.get(0);
    LinkedHashMap<String,String> mods = (LinkedHashMap<String,String>)elModCharge.get(1);
    Hashtable<String,Integer> modToCharge = (Hashtable<String,Integer>)elModCharge.get(2);
    int rowCount = 1;
    Row outRow = sheet.createRow(rowCount);
    Cell label = outRow.createCell(COLUMN_NAME,HSSFCell.CELL_TYPE_STRING);
    label.setCellValue("Name");
    label.setCellStyle(headerStyle);
    label = outRow.createCell(COLUMN_DBS,HSSFCell.CELL_TYPE_STRING);
    label.setCellValue("dbs");
    label.setCellStyle(headerStyle);
    label = outRow.createCell(COLUMN_MOL,HSSFCell.CELL_TYPE_STRING);
    label.setCellValue(HEADER_MOLECULAR_SPECIES_WITH_DOUBLE_BOND_POSITIONS);
    label.setCellStyle(headerStyle);
    int columnNrFirstElement = COLUMN_FIRST_ELEMENT;
    int elementColumns = 0;
    Hashtable<String,Integer> elementColumnLookup = new Hashtable<String,Integer>();
    for (String element : elements) {
      elementColumnLookup.put(element, columnNrFirstElement+elementColumns);
      label = outRow.createCell(columnNrFirstElement+elementColumns,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue(element);
      label.setCellStyle(headerStyle);
      elementColumns++;
    }
    int firstModColumn = columnNrFirstElement+elementColumns;
    int massColumns = 0;
    for (String mod : mods.keySet()) {
      label = outRow.createCell(firstModColumn+massColumns,HSSFCell.CELL_TYPE_STRING);
      String nameHeader = "mass(form["+mods.get(mod)+"] name["+mod+"]";
      if (modToCharge.get(mod)>1)
        nameHeader += " charge="+modToCharge.get(mod);
      nameHeader += ")";
      label.setCellValue(nameHeader);
      label.setCellStyle(headerStyle);   
      massColumns++;
    }
    label = outRow.createCell(firstModColumn+massColumns,HSSFCell.CELL_TYPE_STRING);
    label.setCellValue("tR (min)");
    label.setCellStyle(headerStyle);
    rowCount++;
    
    for (String analyte : analytes) {
      
//      if (labelLookup.containsKey(analyte)) {
//        System.out.println("! "+analyte);
//        for(String labelWhole : labelLookup.get(analyte).keySet()) {
//          System.out.println(labelWhole+":  "+labelLookup.get(analyte).get(labelWhole));
//        }
//      }
      //get out the original results from the labeled samples
      //String molName;
           
      Hashtable<String,QuantVO> quantAnal = quantObjects.get(analyte);
      int modCount = 0;
      outRow = sheet.createRow(rowCount);
      Hashtable<String,Integer> formula = new Hashtable<String,Integer>();
      for (String mod : mods.keySet()) {
        QuantVO quant = quantAnal.get(mod);
        if (modCount==0) {
          label = outRow.createCell(COLUMN_NAME,HSSFCell.CELL_TYPE_STRING);
          label.setCellValue(quant.getAnalyteName());
          label = outRow.createCell(COLUMN_COLON,HSSFCell.CELL_TYPE_STRING);
          label.setCellValue(":");
          label = outRow.createCell(COLUMN_DBS,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(quant.getDbs());
          formula = StaticUtils.categorizeFormula(quant.getAnalyteFormula());
          for (String element : formula.keySet()) {
            label = outRow.createCell(elementColumnLookup.get(element),HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(formula.get(element));
          }
        }
        label = outRow.createCell(firstModColumn+modCount,HSSFCell.CELL_TYPE_NUMERIC);
        label.setCellValue(quant.getAnalyteMass());
        modCount++;
      }
      rowCount++;
      
      if (labeled.containsKey(analyte)) {
        rtsSpeciesAndOmega = new Hashtable<String,Hashtable<String,List<Double>>>();
        for (Integer omegaPos : sameOmegaLabels.keySet()) {
          Hashtable<String,IsotopicLabelVO> labels = sameOmegaLabels.get(omegaPos);
          sumSpeciesHash = new Hashtable<String,Hashtable<Double,DoubleIntegerVO>>();
          molSpeciesHash = new Hashtable<String,Hashtable<String,Hashtable<Double,DoubleIntegerVO>>>();
          rtsOfOneSpecies = new Hashtable<String,Hashtable<String,List<Double>>>(); 
          for(String labelId : labels.keySet()) {
            if (!labeled.get(analyte).containsKey(labelId))
              continue;
            for (AnalyteOmegaInfoVO info :labeled.get(analyte).get(labelId)) {
              //System.out.println(lClass+" "+molecule+": "+info.getLabelDetection()+"     "+label+" ; "+info.getAppliedLabelId());
              for (String exp : analysisModule.getExpNamesInSequence()) {
                areaVO = analysisModule.getResultAreaVO(className,info.getLabelDetection(),exp);
                if (areaVO==null || !areaVO.isMsnEvidenceThere())
                  continue;
                //check if the chemical formula is the one of a labeled species
                elementsLabeled = new Hashtable<String,Integer>(areaVO.getChemicalFormulaElements());
                for (IsotopicLabelVO labelVO : info.getLabels()) {
                  for (String element : labelVO.getLabelElements().keySet()) {
                    nr = 0;
                    if (elementsLabeled.containsKey(element))
                      nr = elementsLabeled.get(element);
                    nr -= labels.get(labelId).getLabelElements().get(element);
                    if (nr==0)
                      elementsLabeled.remove(element);
                    else
                      elementsLabeled.put(element, nr);
                  }
                }
                if (!StaticUtils.isChemicalFormulaTheSame(formula,elementsLabeled))
                  continue;
                //System.out.println("1. "+className+" "+analyte+": "+info.getLabelDetection()+" "+"     "+info.getAppliedLabelId()+" ; "+info.getAppliedLabelId());
                if (areaVO==null || !areaVO.isMsnEvidenceThere())
                  continue;
                if (areaVO.getStrongestChainIdentification()!=null && nrChains>1) {
                  for (String molSpecies : areaVO.getChainInformationTotal().keySet()) {
                    chainsWithoutLabel = StaticUtils.removeLabelsFromChains(molSpecies,info.getLabels());
                    labeledChainsOnly = StaticUtils.getLabeledChainsOnly(molSpecies,info.getLabels());
                    if (labeledChainsOnly==null || labeledChainsOnly.length()==0)
                      continue;
                    encodedOmegaChainInformation = StaticUtils.encodeOmegaPositions(molSpecies,info.getLabels(),analysisModule.getFaHydroxyEncoding(), analysisModule.getLcbHydroxyEncoding());
//                    if (className.equalsIgnoreCase("PC") && encodedOmegaChainInformation.equalsIgnoreCase("16:1(n-7)_18:1"))
//                      System.out.println("chainsWithLabel: "+molSpecies+" - omegaEncoded: "+encodedOmegaChainInformation);
                    if (!molSpeciesHash.containsKey(chainsWithoutLabel))
                      molSpeciesHash.put(chainsWithoutLabel,new Hashtable<String,Hashtable<Double,DoubleIntegerVO>>());
                    if (!molSpeciesHash.get(chainsWithoutLabel).containsKey(encodedOmegaChainInformation))
                      molSpeciesHash.get(chainsWithoutLabel).put(encodedOmegaChainInformation, new Hashtable<Double,DoubleIntegerVO>());
                    //System.out.println("originalRt: "+areaVO.getRtOriginal());
                    rt = Calculator.FormatNumber(info.calculateRtIncludingTheShift(Double.parseDouble(areaVO.getRtOriginal())),5d);
                    //System.out.println("rt-shifted: "+rt);
                    countRts = new DoubleIntegerVO(0,0d);
                    if (molSpeciesHash.get(chainsWithoutLabel).get(encodedOmegaChainInformation).containsKey(rt))
                      countRts = molSpeciesHash.get(chainsWithoutLabel).get(encodedOmegaChainInformation).get(rt);
                    countRts.setKey(countRts.getKey()+1);
                    countRts.setValue(countRts.getValue()+areaVO.getChainInformationTotal().get(molSpecies));
                    molSpeciesHash.get(chainsWithoutLabel).get(encodedOmegaChainInformation).put(rt, countRts);
                  }
                } else if (nrChains==1) {
                  encodedOmegaChainInformation = StaticUtils.encodeOmegaPositions(StaticUtils.encodeLipidCombi(StaticUtils.decodeFAsFromHumanReadableName(info.getLabelDetection().substring(0,info.getLabelDetection().lastIndexOf("_")),analysisModule.getFaHydroxyEncoding(),analysisModule.getLcbHydroxyEncoding(),false,null)),info.getLabels(),analysisModule.getFaHydroxyEncoding(), analysisModule.getLcbHydroxyEncoding());
                  if (!sumSpeciesHash.containsKey(encodedOmegaChainInformation))
                    sumSpeciesHash.put(encodedOmegaChainInformation, new Hashtable<Double,DoubleIntegerVO>());
                  rt = Calculator.FormatNumber(info.calculateRtIncludingTheShift(Double.parseDouble(areaVO.getRtOriginal())),5d);
                  countRts = new DoubleIntegerVO(0,0d);
                  if (sumSpeciesHash.get(encodedOmegaChainInformation).containsKey(rt))
                    countRts = sumSpeciesHash.get(encodedOmegaChainInformation).get(rt);
                  countRts.setKey(countRts.getKey()+1);
                  countRts.setValue(countRts.getValue()+areaVO.getTotalArea(Integer.MAX_VALUE));
                  sumSpeciesHash.get(encodedOmegaChainInformation).put(rt, countRts);
                }
              }
            }
          }
          if (nrChains>1 && molSpeciesHash.size()>0) {
            for (String molSpecies : molSpeciesHash.keySet()) {
              if (molSpeciesHash.get(molSpecies).size()==0)
                continue;
              for (String labeledPart : molSpeciesHash.get(molSpecies).keySet()) {
//                if (className.equalsIgnoreCase("PC") && labeledPart.equalsIgnoreCase("16:1(n-7)_18:1"))
//                  System.out.println("chainsWithLabel: "+molSpecies+" - omegaEncoded: "+labeledPart);
                Vector<Double> meanRts = groupByRtAndCalculateMean(molSpeciesHash.get(molSpecies).get(labeledPart),rtTolerance);
                if (meanRts.size()==0)
                  continue;
                if (!rtsOfOneSpecies.containsKey(molSpecies))
                  rtsOfOneSpecies.put(molSpecies, new Hashtable<String,List<Double>>());
                rtsOfOneSpecies.get(molSpecies).put(labeledPart, meanRts);
              }
            }
          }else if (nrChains==1 && sumSpeciesHash.size()>0) {
            for (String omegaName : sumSpeciesHash.keySet()) {
              Vector<Double> meanRts = groupByRtAndCalculateMean(sumSpeciesHash.get(omegaName),rtTolerance);
              if (meanRts.size()==0)
                continue;
              if (!rtsOfOneSpecies.containsKey(analyte))
                rtsOfOneSpecies.put(analyte, new Hashtable<String,List<Double>>());
                //System.out.println("2. "+className+" ; "+analyte+";"+info.getLabelDetection());
              rtsOfOneSpecies.get(analyte).put(omegaName, meanRts);
            }
          }
          //now group according to the (molecular) species name irrespective of the omega position
          for (String species : rtsOfOneSpecies.keySet()) {
            if (!rtsSpeciesAndOmega.containsKey(species))
              rtsSpeciesAndOmega.put(species, new Hashtable<String,List<Double>>());
            //if (!rtsSpeciesAndOmega.get(species).containsKey(omegaPos))
            //TODO: this might be a good position to do a further grouping of the RTs, even if the things come from a different label
            for (String labeledPart : rtsOfOneSpecies.get(species).keySet()){
              List<Double> rtSorted = rtsOfOneSpecies.get(species).get(labeledPart);
              if (rtsSpeciesAndOmega.get(species).containsKey(labeledPart)) {
                rtSorted.addAll(rtsSpeciesAndOmega.get(species).get(labeledPart));
                rtSorted = groupByRtAndCalculateMean(rtSorted,rtTolerance);                
                Collections.sort(rtSorted);
              }
              rtsSpeciesAndOmega.get(species).put(labeledPart, rtSorted);

            }
//            rtsSpeciesAndOmega.get(species).put(omegaPos, rtsOfOneSpecies.get(species));
          }
          
//          System.out.println("n-"+omegaPos);
//          for (String species : rtsOfOneSpecies.keySet()) {
//            System.out.println("2. "+species+": "+rtsOfOneSpecies.get(species));
//          }
          
        }
        
          
        Vector<String> sortedSpecies =  StaticUtils.sortChainCombinations(rtsSpeciesAndOmega.keySet());
        for (String species : sortedSpecies) {
          Hashtable<String,List<Double>> rtsOfSpecies = rtsSpeciesAndOmega.get(species);
          List<Double> sortedRts = new ArrayList<Double>();
          Hashtable<Double,Double> used = new Hashtable<Double,Double>();
          for (List<Double> rts : rtsOfSpecies.values()) {
            for (Double oneRt : rts) {
              if (used.containsKey(oneRt))
                continue;
              sortedRts.add(oneRt);
              used.put(oneRt, oneRt);
            }
          }
          Collections.sort(sortedRts);
          
          
          // this checks if a hit with less specified omega positions matches a hit with more
          Double lastRt = new Double (0);
          Hashtable<String,Hashtable<Double,Double>> leaveItOut = new Hashtable<String,Hashtable<Double,Double>>();
          for (Double oneRt : sortedRts) {
            if ((oneRt-lastRt)<labelTol) {
              for (String omegaEncoded1: rtsOfSpecies.keySet()) {
                if (!rtsOfSpecies.get(omegaEncoded1).contains(lastRt))
                  continue;
                for (String omegaEncoded2: rtsOfSpecies.keySet()) {
                  if (!rtsOfSpecies.get(omegaEncoded2).contains(oneRt))
                    continue;
                  String obsoleteOne = findAnObsoleteIdentification(omegaEncoded1,omegaEncoded2,analysisModule.getFaHydroxyEncoding(),analysisModule.getLcbHydroxyEncoding());
                  if (obsoleteOne==null)
                    continue;
                  if (!leaveItOut.containsKey(obsoleteOne))
                    leaveItOut.put(obsoleteOne, new Hashtable<Double,Double>());
                  if (obsoleteOne.equalsIgnoreCase(omegaEncoded1))
                    leaveItOut.get(obsoleteOne).put(lastRt, lastRt);
                  else
                    leaveItOut.get(obsoleteOne).put(oneRt, oneRt);
                }
              }
            }
            lastRt = oneRt;
          }
          
          for (Double oneRt : sortedRts) {
            for (String omegaEncoded: rtsOfSpecies.keySet()) {
              if (!rtsOfSpecies.get(omegaEncoded).contains(oneRt) || (leaveItOut.containsKey(omegaEncoded) && leaveItOut.get(omegaEncoded).containsKey(oneRt)))
                continue;
              modCount = 0;
              outRow = sheet.createRow(rowCount);
              for (String mod : mods.keySet()) {
                QuantVO quant = quantAnal.get(mod);
                if (modCount==0) {
                  label = outRow.createCell(COLUMN_NAME,HSSFCell.CELL_TYPE_STRING);
                  label.setCellValue(quant.getAnalyteName());
                  label = outRow.createCell(COLUMN_COLON,HSSFCell.CELL_TYPE_STRING);
                  label.setCellValue(":");
                  label = outRow.createCell(COLUMN_DBS,HSSFCell.CELL_TYPE_NUMERIC);
                  label.setCellValue(quant.getDbs());
                  label = outRow.createCell(COLUMN_MOL,HSSFCell.CELL_TYPE_STRING);
                  label.setCellValue(omegaEncoded);
                  for (String element : formula.keySet()) {
                    label = outRow.createCell(elementColumnLookup.get(element),HSSFCell.CELL_TYPE_NUMERIC);
                    label.setCellValue(formula.get(element));
                  }
                }
                label = outRow.createCell(firstModColumn+modCount,HSSFCell.CELL_TYPE_NUMERIC);
                label.setCellValue(quant.getAnalyteMass());
                label = outRow.createCell(firstModColumn+massColumns,HSSFCell.CELL_TYPE_NUMERIC);
                label.setCellValue(Calculator.FormatNumber(oneRt.doubleValue(), 2d));

                modCount++;
              }
              rowCount++;

            }
          }
        }
      }
    }
  }
  
  
  /**
   * extracts from the original mass list the chemical elements used, the modifications applied and the charge of those modifications
   * @param quantObjects the information form the original Excel mass list; first key: analyte name; second key modification name; value: information about the analyte entry
   * @return vector containing three objects: first object: a vector containing the elements in the order of the Hill notation; second object: a LinkedHashMap containing the modifications as key and the modification formula as value; third object: a hash table; key: modification name; value: charge of the modification 
   * @throws ChemicalFormulaException thrown when there is something wrong with the chemical formula
   */
  private Vector<Object> getAvailableElementsAndModificationsPlusCharge(Hashtable<String,Hashtable<String,QuantVO>> quantObjects) throws ChemicalFormulaException{
    Vector<Object> result = new Vector<Object>();
    Set<String> elements = new HashSet<String>();
    LinkedHashMap<String,String> modifications = new LinkedHashMap<String,String>();
    Hashtable<String,Integer> modToCharge = new Hashtable<String,Integer>();
    for (Hashtable<String,QuantVO> quantAnal : quantObjects.values()) {
      for (QuantVO quant : quantAnal.values()) {
        for (String element : StaticUtils.categorizeFormula(quant.getAnalyteFormula()).keySet()) {
          if (!elements.contains(element))
            elements.add(element);
          if (!modifications.containsKey(quant.getModName())) {
            modifications.put(quant.getModName(),StaticUtils.getFormulaInHillNotation(StaticUtils.categorizeFormula(quant.getModFormula()),false));
            modToCharge.put(quant.getModName(), quant.getCharge());
          }
        }
      }
    }
    Vector<String> sorted = new Vector<String>();
    if (elements.contains("C"))
      sorted.add("C");
    if (elements.contains("H"))
      sorted.add("H");
    List<String> otherThanCH = new ArrayList<String>();
    for (String element : elements) {
      if (!element.equalsIgnoreCase("C") && !element.equalsIgnoreCase("H"))
        otherThanCH.add(element);
    }
    Collections.sort(otherThanCH);
    for (String element : otherThanCH)
      sorted.add(element);
    result.add(sorted);
    result.add(modifications);
    result.add(modToCharge);
    return result;
  }
  
  
  private static XSSFCellStyle getHeaderStyle(XSSFWorkbook wb){
    XSSFCellStyle arial12style = wb.createCellStyle();
    XSSFFont arial12font = wb.createFont();
    arial12font.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(HSSFCellStyle.ALIGN_CENTER);
    return arial12style;
  }
  
  
  /**
   * calculates the standard deviation of the analytes with one double bond
   * @param acceptedMolecules the accepted labeled molecules
   * @param analysisModule the results read from the LDA files
   * @return the standard deviation
   */
  private double getMeanStandardDeviationOfSingleDoubleBonds(Hashtable<String,Hashtable<String,Hashtable<String,Vector<AnalyteOmegaInfoVO>>>> acceptedMolecules, ComparativeAnalysis analysisModule) {   
    double sumSquaredDevs = 0d;
    int count = 0;
    int nrChains;
    ResultAreaVO resultVO;
    //first key: expName; second key: RT
    Hashtable<String,Hashtable<String,ResultAreaVO>> sumSpeciesHash;
    //first key: chain idents; second key: expName; third key: RT
    Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>> molSpeciesHash;
    for (String lClass : acceptedMolecules.keySet()) {
      nrChains = analysisModule.getNrOfChainsOfClass().get(lClass);
      for (String molecule : acceptedMolecules.get(lClass).keySet()) {
        for (String label : acceptedMolecules.get(lClass).get(molecule).keySet()) {
          sumSpeciesHash = new Hashtable<String,Hashtable<String,ResultAreaVO>>();
          molSpeciesHash = new Hashtable<String,Hashtable<String,Hashtable<String,ResultAreaVO>>>();
          for (AnalyteOmegaInfoVO info : acceptedMolecules.get(lClass).get(molecule).get(label)) {
            //System.out.println(lClass+" "+molecule+": "+info.getLabelDetection()+"     "+label+" ; "+info.getAppliedLabelId());
            if (!info.isSingelDoubleBond())
              continue;
            for (String exp : analysisModule.getExpNamesInSequence()) {
              resultVO = analysisModule.getResultAreaVO(lClass,info.getLabelDetection(),exp);
              if (resultVO==null || !resultVO.isMsnEvidenceThere())
                continue;
              if (resultVO.getStrongestChainIdentification()!=null && nrChains>1) {
                for (String molSpecies : resultVO.getChainInformationTotal().keySet()) {
                  if (!molSpeciesHash.containsKey(molSpecies))
                    molSpeciesHash.put(molSpecies,new Hashtable<String,Hashtable<String,ResultAreaVO>>());
                  if (!molSpeciesHash.get(molSpecies).containsKey(exp))
                    molSpeciesHash.get(molSpecies).put(exp, new Hashtable<String,ResultAreaVO>());
                  molSpeciesHash.get(molSpecies).get(exp).put(resultVO.getRtOriginal(), resultVO);
                }
              }else if (nrChains==1) {
                if (!sumSpeciesHash.containsKey(exp))
                  sumSpeciesHash.put(exp, new Hashtable<String,ResultAreaVO>());
                sumSpeciesHash.get(exp).put(resultVO.getRtOriginal(), resultVO);
              }
            }
          }
          
          if (molSpeciesHash.size()>0) {
            for (String molSpecies : molSpeciesHash.keySet()) {
              if (!moreThanOneHit(molSpeciesHash.get(molSpecies)))
                continue;
              //System.out.println(molSpecies);
              Vector<Double> squareDevs = calculateSquareDeviationValuesOfWeightedMean(molSpeciesHash.get(molSpecies),molSpecies);
              count += squareDevs.size();
              for (Double squareDev : squareDevs)
                sumSquaredDevs += squareDev;
            }
          } else {
            if (moreThanOneHit(sumSpeciesHash)) {
              //System.out.println("SUM!!!");
              Vector<Double> squareDevs = calculateSquareDeviationValuesOfWeightedMean(sumSpeciesHash,null);
              count += squareDevs.size();
              for (Double squareDev : squareDevs)
                sumSquaredDevs += squareDev;

            }
          }
        }
        //System.out.println("Accepted: "+lClass+" "+molecule);
      }
    }
    return Math.sqrt(sumSquaredDevs)/((double)(count-1));
  }


  /**
   * checks whether more than one hit has been detected for this analyte
   * @param results the results to check; first key: exp name; second key: RT; value the VO containing information on the detection
   * @return true when there is more than one hit present
   */
  private boolean moreThanOneHit(Hashtable<String,Hashtable<String,ResultAreaVO>> results) {
    int count = 0;
    for (Hashtable<String,ResultAreaVO> hits : results.values()) {
      count+=hits.size();
      if (count>1)
        return true;
    }
    return false;
  }
  

  /**
   * first, calculates a mean weighted by the areas of detection; second, calculates the squared deviations from this weighted mean
   * @param results the results to check; first key: exp name; second key: RT; value the VO containing information on the detection
   * @param molSpecies the molecular species
   * @return the squared deviation values from the weighted mean to calculate a standard deviation
   */
  private Vector<Double> calculateSquareDeviationValuesOfWeightedMean(Hashtable<String,Hashtable<String,ResultAreaVO>> results, String molSpecies){
    //first, calculated the weighted mean
    double sum = 0d;
    double sumWeighting = 0d;
    for (Hashtable<String,ResultAreaVO> hits : results.values()) {
      for (String rtString : hits.keySet()) {
        double rt = Double.valueOf(rtString);
        double area = molSpecies != null ? hits.get(rtString).getChainInformationTotal().get(molSpecies) : hits.get(rtString).getTotalArea(hits.get(rtString).getMaxIsotope());
        sum += rt*area;
        sumWeighting += area;
      }
    }
    double weightedMean = sum/sumWeighting;
    //System.out.println("weightedMean: "+weightedMean);
    Vector<Double> devs = new Vector<Double>();
    for (Hashtable<String,ResultAreaVO> hits : results.values()) {
      for (String rtString : hits.keySet()) {
        double rt = Double.valueOf(rtString);
        devs.add(Math.pow(rt-weightedMean, 2d));
      }
    }
    return devs;
  }
  
  /**
   * groups hits of different retention times; the rtTolerance value dictates the maximum tolerance; when a group is defined the mean value of the retention times is calculated and used as reference value;
   * when there are several options, according to the settings, this method returns all values (OMEGA_PEAK_RTDIFF) or the strongest identification only (OMEGA_PEAK_STRONGEST)
   * @param rtsIn a hash table containing the retention times and information about this hit; key: retention time; value: DoubleIntegerVO where the key is the number of occurrences of exactly this RT, and the value is the total area
   * @param rtTolerance the maximum RT tolerance for assigning a hit to a group
   * @return the retention time(s) of the (strongest) group(s)
   */
  private Vector<Double> groupByRtAndCalculateMean(Hashtable<Double,DoubleIntegerVO> rtsIn, double rtTolerance){
    Vector<Double> meanRts = groupByRtAndCalculateMean(rtsIn.keySet(), rtTolerance);
    if (Settings.getOmegaExportPeakSelection()==Settings.OMEGA_PEAK_RTDIFF)
      return meanRts;
    else if (Settings.getOmegaExportPeakSelection()==Settings.OMEGA_PEAK_STRONGEST) {
      Hashtable<Double,Double> areas = new Hashtable<Double,Double>();
      for (Double rt : meanRts)
        areas.put(rt, 0d);
      for (Double other : rtsIn.keySet()) {
        for (Double rt : areas.keySet()) {
          if (Math.abs(rt-other)<rtTolerance)
            areas.put(rt, (areas.get(rt)+rtsIn.get(other).getValue()));
        }
      }
      double highest = 0d;
      Double rtHighest = 0d;
      for (Double rt : areas.keySet()) {
        if (areas.get(rt)>highest) {
          highest = areas.get(rt);
          rtHighest = rt;
        }
      }
      meanRts = new Vector<Double>();
      meanRts.add(rtHighest);
    }
    return meanRts;
  }
  
  /**
   * groups hits of different retention times; the rtTolerance value dictates the maximum tolerance; when a group is defined the mean value of the retention times is calculated and used as reference value;
   * @param rtsIn a list of retention times
   * @param rtTolerance he maximum RT tolerance for assigning a hit to a group
   * @return the retention time(s) of the (strongest) group(s)
   */
  private Vector<Double> groupByRtAndCalculateMean(Collection<Double> rtsIn, double rtTolerance){
    List<Double> rtSorted = new ArrayList<Double>(rtsIn);
    Vector<Vector<Double>> rtGroups = new Vector<Vector<Double>>();
    Collections.sort(rtSorted);
    for (Double oneRt : rtSorted) {
      if (rtGroups.size()==0 || (oneRt-rtGroups.lastElement().lastElement())>rtTolerance)
        rtGroups.add(new Vector<Double>());
      rtGroups.lastElement().add(oneRt);
    }
    Vector<Double> meanRts = new Vector<Double>();
    if (rtGroups.size()==0)
      return meanRts;
    for (Vector<Double> rts : rtGroups) {
      double sum = 0d;
      for (double oneRt : rts)
        sum += oneRt;
      meanRts.add(sum/((double)rts.size()));
    }
    return meanRts;
  }
  
  /**
   * checks whether a hit coincides with a hit where more omega positions are determined
   * @param one the first molecular species identification
   * @param two the second molecular species identification
   * @param faHydroxyEncoding the OH encodings of the FA moiety
   * @param lcbHydroxyEncoding the OH encodings of the LCB moiety
   * @return the id if the obsolete hit, or null if none of the hits are obsolete
   * @throws LipidCombinameEncodingException when there is something wrong with the lipid name encoding
   */
  private String findAnObsoleteIdentification(String one, String two, HydroxyEncoding faHydroxyEncoding, HydroxyEncoding lcbHydroxyEncoding) throws LipidCombinameEncodingException {
    Vector<FattyAcidVO> first = StaticUtils.decodeFAsFromHumanReadableName(one, faHydroxyEncoding, lcbHydroxyEncoding, false, null);
    Vector<FattyAcidVO> second = StaticUtils.decodeFAsFromHumanReadableName(two, faHydroxyEncoding, lcbHydroxyEncoding, false, null);
    if (first.size()!=second.size())
      return null;
    int omegaPos = -1;
    int nrNFirst = 0;
    int nrNSecond = 0;
    for (FattyAcidVO fa : first) {
      if (fa.getOmegaPosition()>-1) {
        nrNFirst++;
        if (omegaPos==-1)
          omegaPos=fa.getOmegaPosition();
        else if (omegaPos!=fa.getOmegaPosition())
          return null;
      }
    }
    for (FattyAcidVO fa : second) {
      if (fa.getOmegaPosition()>-1) {
        nrNSecond++;
        if (omegaPos==-1)
          omegaPos=fa.getOmegaPosition();
        else if (omegaPos!=fa.getOmegaPosition())
          return null;
      }
    }
    if (nrNFirst==nrNSecond)
      return null;
    else if (nrNFirst>nrNSecond)
      return two;
    else
      return one;
  }
  
}
