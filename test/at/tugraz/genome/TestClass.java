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

package at.tugraz.genome;

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.LineNumberReader;
import java.io.OutputStream;
import java.math.BigDecimal;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TimerTask;
import java.util.Vector;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.imageio.ImageIO;
import javax.swing.JApplet;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JProgressBar;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFFont;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.util.CellRangeAddress;
import org.apache.poi.xssf.usermodel.XSSFCell;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFColor;
import org.apache.poi.xssf.usermodel.XSSFFont;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
//import org.jfree.chart.ChartFactory;
//import org.jfree.chart.JFreeChart;
//import org.jfree.chart.plot.PlotOrientation;
//import org.jfree.chart.plot.Zoomable;
//import org.jfree.data.xy.XYSeries;
//import org.jfree.data.xy.XYSeriesCollection;
//import org.jfree.chart.ChartFactory;
//import org.jfree.chart.JFreeChart;
//import org.jfree.chart.plot.PlotOrientation;
//import org.jfree.chart.plot.XYPlot;
//import org.jfree.chart.renderer.xy.XYBarRenderer;
//import org.jfree.data.xy.XYSeries;
//import org.jfree.data.xy.XYSeriesCollection;
import org.xml.sax.SAXException;













import de.isas.mztab2.io.MzTabFileParser;
import de.isas.lipidomics.mztab2.validation.MzTabValidator;
import uk.ac.ebi.pride.jmztab2.utils.errors.MZTabError;
import uk.ac.ebi.pride.jmztab2.utils.errors.MZTabErrorList;
import uk.ac.ebi.pride.jmztab2.utils.errors.MZTabErrorType;
////import at.tugraz.genome.IndependentTwoSamplesTTest;
import JSci.maths.statistics.ChiSqrDistribution;
import JSci.maths.statistics.TDistribution;
import at.tugraz.genome.dbutilities.SimpleValueObject;
import at.tugraz.genome.exception.LipidBLASTException;
import at.tugraz.genome.exception.MSDialException;
import at.tugraz.genome.lda.BatchQuantThread;
import at.tugraz.genome.lda.LDAResultReader;
import at.tugraz.genome.lda.LipidDataAnalyzer;
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.MzxmlToChromThread;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.QuantificationThread;
import at.tugraz.genome.lda.SingleQuantThread;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.alex123.TargetlistDirParser;
import at.tugraz.genome.lda.alex123.TargetlistParser;
import at.tugraz.genome.lda.alex123.vos.TargetlistEntry;
import at.tugraz.genome.lda.exception.AlexTargetlistParserException;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.HydroxylationEncodingException;
import at.tugraz.genome.lda.exception.LMException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.RulesException;
import at.tugraz.genome.lda.msn.FragmentCalculator;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.MSnAnalyzer;
import at.tugraz.genome.lda.msn.PostQuantificationProcessor;
import at.tugraz.genome.lda.msn.RulesContainer;
import at.tugraz.genome.lda.msn.parser.FALibParser;
import at.tugraz.genome.lda.msn.parser.FragRuleParser;
import at.tugraz.genome.lda.msn.parser.LCBLibParser;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.msn.vos.IntensityChainVO;
import at.tugraz.genome.lda.msn.vos.IntensityPositionVO;
import at.tugraz.genome.lda.msn.vos.IntensityRuleVO;
import at.tugraz.genome.lda.msn.vos.MSnDebugVO;
import at.tugraz.genome.lda.parser.MzXMLMergerForWaters;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.LipidomicsChromatogram;
import at.tugraz.genome.lda.quantification.LipidomicsDefines;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.quantification.SavGolJNI;
import at.tugraz.genome.lda.swing.AbsoluteQuantSettingsPanel;
import at.tugraz.genome.lda.swing.BarChartPainter;
import at.tugraz.genome.lda.swing.BatchQuantificationTable;
import at.tugraz.genome.lda.swing.BatchQuantificationTableModel;
import at.tugraz.genome.lda.swing.ExportPanel;
import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.LMQuadraticTwoVariables;
import at.tugraz.genome.lda.utils.LevenbergMarquardtOptimizer;
import at.tugraz.genome.lda.utils.RangeInteger;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.lda.vos.RawQuantificationPairVO;
import at.tugraz.genome.lda.xml.AbsoluteQuantSettingsWholeReader;
import at.tugraz.genome.lda.xml.AddScan;
import at.tugraz.genome.lda.xml.MzXmlReader;
import at.tugraz.genome.lda.xml.RawToChromTranslator;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.parser.spectrummill.vos.AaFormulaVO;
import at.tugraz.genome.maspectras.parser.spectrummill.vos.ChemicalFormulaVO;
import at.tugraz.genome.maspectras.parser.spectrummill.vos.SmChemicalElementVO;
import at.tugraz.genome.maspectras.parser.spectrummill.vos.SmIsotopeVO;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.maspectras.quantification.CgChromatogram;
import at.tugraz.genome.maspectras.quantification.CgDefines;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgIAddScan;
import at.tugraz.genome.maspectras.quantification.CgMzXmlReader;
import at.tugraz.genome.maspectras.quantification.CgParameterSet;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.CgReader;
import at.tugraz.genome.maspectras.quantification.CgScan;
import at.tugraz.genome.maspectras.quantification.CgScanHeader;
import at.tugraz.genome.maspectras.quantification.ChromatogramReader;
import at.tugraz.genome.maspectras.quantification.Probe3D;
import at.tugraz.genome.maspectras.quantification.RawToChromatogramTranslator;
import at.tugraz.genome.maspectras.utils.Calculator;
import at.tugraz.genome.maspectras.utils.StringUtils;
import at.tugraz.genome.parsers.LipidBLASTParser;
import at.tugraz.genome.parsers.MSDialTxtParser;
import at.tugraz.genome.parsers.MSFinderStructureParser;
import at.tugraz.genome.util.FloatMatrix;
import at.tugraz.genome.util.index.IndexFileException;
import at.tugraz.genome.vos.BrainSpecies;
import at.tugraz.genome.vos.FoundBiologicalSpecies;
import at.tugraz.genome.vos.LdaDialStandardsEvidence;
import at.tugraz.genome.vos.LdaLBLASTCompareVO;
import at.tugraz.genome.vos.LdaLbStandardsEvidence;
import at.tugraz.genome.vos.LipidBLASTDetectionVO;
import at.tugraz.genome.vos.LipidBLASTIdentificationVO;
import at.tugraz.genome.vos.LipidClassInfoVO;
import at.tugraz.genome.vos.MSDialEntry;
import at.tugraz.genome.vos.MSFinderEntry;
import at.tugraz.genome.vos.MSFinderHitVO;
import at.tugraz.genome.vos.Mix1Standards;
import at.tugraz.genome.vos.Mix2Standards;
import at.tugraz.genome.vos.ReferenceInfoVO;
import at.tugraz.genome.voutils.GeneralComparator;
//import at.tugraz.genome.thermo.IXRawfile;
//import at.tugraz.genome.thermo.IXRawfile3;




//import javax.imageio.ImageIO;
//
//import org.jfree.chart.ChartFactory;
//import org.jfree.chart.JFreeChart;
//import org.jfree.chart.axis.LogarithmicAxis;
//import org.jfree.chart.labels.XYToolTipGenerator;
//import org.jfree.chart.plot.PlotOrientation;
//import org.jfree.chart.plot.XYPlot;
//import org.jfree.chart.renderer.xy.XYBarRenderer;
//import org.jfree.chart.renderer.xy.XYItemRenderer;
//import org.jfree.data.xy.XYDataItem;
//import org.jfree.data.xy.XYDataset;
//import org.jfree.data.xy.XYSeries;
//import org.jfree.data.xy.XYSeriesCollection;



























import com.sun.j3d.utils.applet.MainFrame;
import com.sun.org.apache.xerces.internal.util.URI;

import de.isas.mztab2.model.Metadata;
import de.isas.mztab2.model.MzTab;
import de.isas.mztab2.model.ValidationMessage;

//import com.sun.jna.Native;
//import com.sun.jna.NativeLibrary;
//import com4j.Holder;
//import com4j.Variant;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class TestClass extends JApplet implements AddScan
{
  public static void main(String[] args)
  {
    // MainFrame frame = new MainFrame(new TestClass(), 1024, 1024);
    new TestClass();
  }
  
  private final static int NOT_FOUND = 0;
  private final static int MS1_FOUND = 1;
  private final static int FA_FOUND = 2;
  private final static int FA_CORRECT = 3;
  private final static int POS_CORRECT = 4;
  
  private final static double RT_TOL = 0.6d;
  
  private final static double MSDIAL_MZ_TOL = 0.006d;
  private final static float MSDIAL_RT_TOL = 0.25f;

  public TestClass()
  {
    //int charge = 2;
    //double number = 0.000000000000001/(float)charge;
    //System.out.println(number);
    // this.testPaintBarChart();
    // String[] names = ImageIO.getWriterFormatNames();
    // for (String name:names)
    // System.out.println(name);
    //this.extractSameMasses();
    // this.testExportPanel();
    //this.translateTAGListToMzMineFormat();
    //this.shortenMSList();
    //this.calculateTheoreticalMass();
    //calculateByKnownMassesAndCAtoms();
    //this.calculateIntensityDistribution();
    //this.justMzValues();
    //this.ttest();
    //this.ttestFromExcel();
    //this.msFileReaderAccess();
    //paintSpectrumFromXls();
    //this.plus2IsotopicRelationXls();
    //this.calculateIntensityDistribution();
    //this.readQuantExcel();
    //this.translateMarleneFileToQuant();
    //this.mergeExcelResults();
    //this.calcPotentialFormulasPerMass();
    //this.translateSabineFileToQuant();
    //this.translateMzMineFilesToExcel();
    //this.mergeNamedAndUntargetedExcel();
    //this.detectMsn();
    //this.readMSnIdentification();
   //this.sumVarianceErrorPropagated();
//    try {
//      this.testTabFile();
//    }
//    catch (Exception e) {
//      // TODO Auto-generated catch block
//      e.printStackTrace();
//    }
//    createN15MassList();
    //this.calculateErrorPropagatedStdev();
    //calculateRetentionTime();
    //paintRtDependencies();
    //leastSquareFittingFor2Variables();
    
    //findPCAMostContributingComponents();
    //paintPCAOfMostProminentComponents();
    
    //makeCommaSeparatedFiles();
    //doPostQuantProcessing();
    //mergeLipidClassesAndCalcPercentage();
    //this.generateAlkylLinkedMassList();
    //this.detectPlasmalogens();
    //this.alkylAlkenylCombinations();
    //this.searchMSnSpectra();
    //this.mergeMzXMLFiles();
    //this.convertWiffFile();
    //this.evaluateExperiment3();
    //this.writeLDAResultsToCompareExcel();
    //this.compareLDABLASTNaturalProbesNegative();
    //this.compareLDABLASTNaturalProbesPositive();
    //this.compareLDABLASTControlledPositiveProbes();
    //this.compareLDABLASTControlledNegativeProbes();
    //this.testAfterStructuralIdentification();
    //this.quantifyPeak();
    //this.checkForMassDeviation();
    //this.parseRule();
    //this.countMS2();
    //this.mergeTGRessults();
    //this.detectLBNotDetected();
    //this.tryXmlStax();
    //this.generateDetailsBiologicalExperiment();
    //this.generateDetailsExperiment1();
    //this.generateDetailsExperiment2();
    //this.generateDetailsExperiment3();
    //this.calcMeanAndStdev();
    //this.generateDetectedSpeciesListForTG4000QTRAP();
    //this.crossPlatformComparison();
    //this.positionalIsomersSmallerRange();
    //compareDeviationFromTheoreticalValueBasedOnIntensity();
    //this.readAlex123File();
    //this.startQuantitationWithAlex123File();
    //this.readIndexFile();
    //this.readRttFile();
    //this.detectMsnByAlex123();
    //this.mergeIdx2();
    //this.faSummaryConverter();
    //this.mzTabValidation();
    //this.evaluatePrmData();
    //this.batchQuantByCommandLine();
    //this.groupAlexResultsByRt();
    //generateSphingolipidMassList();
    //parseLCBList();
    //replicateCudaBug();
    detectMSnWithLCB();
  }

  private void testExportPanel()
  {
    this.add(new ExportPanel(Color.WHITE, Color.BLACK, null));
  }

  private void testPaintBarChart()
  {
    // BarChartPainter painter = new BarChartPainter("54:6");
    // AbsoluteQuantSettingsPanel panel = new AbsoluteQuantSettingsPanel();
    // this.add(panel);
    // panel.showSettingsPanel();

    try {
      AbsoluteQuantSettingsWholeReader reader = new AbsoluteQuantSettingsWholeReader(
          "X:\\ForJuergen\\ERA000159.sample2.xml", null);
    }
    catch (SAXException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    // try {
    // LineNumberReader reader = new LineNumberReader(new
        // FileReader("X:\\ForJuergen\\ERA000159.sample.xml"));
    // String line;
    // String output = "";
    // while ((line = reader.readLine()) != null) {
    // char myChar = '\u0010';
    // output += line.replaceAll(String.valueOf(myChar), " ");
    // System.out.println(reader.getLineNumber());
    // }
    // reader.close();
    // BufferedOutputStream out = new BufferedOutputStream(new
        // FileOutputStream("X:\\ForJuergen\\ERA000159.sample2.xml"));
    // out.write(output.getBytes());
    // out.close();
    // }
    // catch (FileNotFoundException e) {
    // // TODO Auto-generated catch block
    // e.printStackTrace();
    // } catch (IOException e) {
    // // TODO Auto-generated catch block
    // e.printStackTrace();
    // }
  }

  private void extractSameMasses()
  {
    try {
      Hashtable<Integer,Hashtable<String,Float>> sameMasses = new Hashtable<Integer,Hashtable<String,Float>>();
      InputStream myxls = new FileInputStream("E:\\lipidomics\\20131029\\UNGERADE_pos.xls");
      HSSFWorkbook workbook     = new HSSFWorkbook(myxls);
      float toWholeNumberValue = 100f;
      for (int sheetNumber = 0; sheetNumber!=workbook.getNumberOfSheets(); sheetNumber++){
        HSSFSheet sheet = workbook.getSheetAt(sheetNumber); 
        int sideChainColumn = -1;
        int doubleBondColumn = -1;

        Hashtable<Integer,String> massOfInterestColumns = new Hashtable<Integer,String>();
        Hashtable<String,Hashtable<String,Integer>> adductComposition = new Hashtable<String,Hashtable<String,Integer>>();
        int retTimeColumn = -1;
        boolean foundColumns = false;
        ElementConfigParser aaParser = Settings.getElementParser();
        float fixedStartTime = 0;
        float fixedEndTime = Float.MAX_VALUE;
        Hashtable<Integer,String> elementColumns = new  Hashtable<Integer,String>();
        for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
          HSSFRow row = sheet.getRow(rowCount);
          String sideChain = "";
          int doubleBonds = -1;
          
          Hashtable<String,Integer> elementalComposition = new Hashtable<String,Integer>();
          Hashtable <String,Double> massesOfInterest = new Hashtable <String,Double>();
          float retTime = -1;
          Hashtable<Integer,String> possibleElementColumns = new  Hashtable<Integer,String>();
          for (int i=0; row!=null && i!=(row.getLastCellNum()+1);i++){
            HSSFCell cell = row.getCell(i);
            String contents = "";
            Double numeric = null;
            int cellType = -1;
            if (cell!=null) cellType = cell.getCellType();
            if (cellType==HSSFCell.CELL_TYPE_STRING){
              contents = cell.getStringCellValue();
              try{ numeric = new Double(contents);}catch(NumberFormatException nfx){};
            }else if (cellType==HSSFCell.CELL_TYPE_NUMERIC || cellType==HSSFCell.CELL_TYPE_FORMULA){
             numeric = cell.getNumericCellValue();
             contents = String.valueOf(numeric);
            }
            //String contents = sheet.getCell(i,rowCount).getContents();
            if (contents!=null)
              contents = contents.trim();
            if (!foundColumns){
              if (contents.equalsIgnoreCase("Seitenkette")||contents.equalsIgnoreCase("Name")){
                sideChainColumn = i;
              } else if (contents.equalsIgnoreCase("dbs")||contents.equalsIgnoreCase("dbs_TAG")){
                doubleBondColumn = i;
              } 

              else if (contents.startsWith("mass")&&contents.contains("(")&&contents.contains(")")){
                String[] formulaAndName = QuantificationThread.extractFormulaAndAdductName(contents);
                adductComposition.put(formulaAndName[1],StaticUtils.categorizeFormula(formulaAndName[0]));
                massOfInterestColumns.put(i,formulaAndName[1]);
              }
              else if (contents.equalsIgnoreCase("tR (min)")){
                retTimeColumn = i;
              }
              else if (contents.startsWith("Start-RT:")){
                try{
                  fixedStartTime = Float.parseFloat(contents.substring("Start-RT:".length()).trim().replaceAll(",", "."));
                }catch(NumberFormatException nfx){nfx.printStackTrace();};
              }
              else if (contents.startsWith("Stop-RT:")){
                try{
                  fixedEndTime = Float.parseFloat(contents.substring("Stop-RT:".length()).trim().replaceAll(",", "."));
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
                  if (ok)possibleElementColumns.put(i, contents.trim());
                }
              }
            }else{            
              if (i==sideChainColumn&&contents!=null&contents.length()>0){
//                if (numeric!=null)
//                  sideChain = cell.getStringCellValue();
//                else
                  sideChain = contents;
              }
              if (i==doubleBondColumn&&contents!=null&&contents.length()>0){
                doubleBonds = numeric.intValue();
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
//                if (contents.contains(","))
//                  contents = contents.replaceAll(",", ".");
                //massOfInterest = Float.parseFloat(contents);
                massesOfInterest.put(massOfInterestColumns.get(i), numeric);
              }
              if (i==retTimeColumn&&contents!=null&contents.length()>0){
//                if (contents.contains(","))
//                  contents = contents.replaceAll(",", ".");
                retTime = numeric.floatValue();
              }
            }
          }

          if (sideChainColumn != -1 && possibleElementColumns.size()>0 && 
              massOfInterestColumns.size()>0 && retTimeColumn != -1&&
              !foundColumns){
            foundColumns = true;
            elementColumns = new Hashtable<Integer,String>(possibleElementColumns);
          }
          for (String modName : massesOfInterest.keySet()){
            double massOfInterest = massesOfInterest.get(modName);
            String doubleString = String.valueOf(Calculator.roundFloat((float)massOfInterest * toWholeNumberValue, 0, BigDecimal.ROUND_DOWN));
            Integer intValue = new Integer(doubleString.substring(0,doubleString.indexOf(".")));
            Hashtable<String,Float> hash = new Hashtable<String,Float>();
            if (sameMasses.containsKey(intValue))
              hash = sameMasses.get(intValue);
            hash.put(sheet.getSheetName() + "_" + StaticUtils.generateLipidNameString(sideChain, doubleBonds)+"_"+modName, (float)massOfInterest);
            sameMasses.put(intValue, hash);
          }
        }
      }
      myxls.close();
      List<Integer> keyList = new ArrayList<Integer>(sameMasses.keySet());
      for (Integer massKey : keyList){
        Vector<Float> allComparableMasses = new Vector<Float>();
        Vector<String> analytesSorted = new  Vector<String>();
        Hashtable<String,Float> sameMass = sameMasses.get(massKey);
        for (String analyte : sameMass.keySet()){
          analytesSorted.add(analyte);
          allComparableMasses.add(sameMass.get(analyte));
        }
        // in order to avoid doubling of entries it is sufficient to check the group that is above this one
        if (sameMasses.containsKey((massKey+1))){
          sameMass = sameMasses.get(massKey+1);
          for (String analyte : sameMass.keySet()){
            analytesSorted.add(analyte);
            allComparableMasses.add(sameMass.get(analyte));
          }
        }
        if (allComparableMasses.size()>1){
          for (int i=0; i!=allComparableMasses.size(); i++){
            Hashtable<String,Float> identifiedSame = new Hashtable<String,Float>();
            identifiedSame.put(analytesSorted.get(i), allComparableMasses.get(i));
            for (int j=i+1; j!=allComparableMasses.size(); j++){
              if ((allComparableMasses.get(i)-(1f/toWholeNumberValue))<=allComparableMasses.get(j)&&allComparableMasses.get(j)<=(allComparableMasses.get(i)+(1f/toWholeNumberValue))){
                identifiedSame.put(analytesSorted.get(j), allComparableMasses.get(j));
              }
            }
            if (identifiedSame.size()>1){
              for (String key:identifiedSame.keySet()) {
                System.out.println(key + "\t" + identifiedSame.get(key));
              }
              System.out.println("--------------------------");
            }
          }
        }  
      }
      
//      for (Hashtable<String,Float> sameMass:sameMasses.values()) {
//        if (sameMass.size() > 1) {
//          for (String key:sameMass.keySet()) {
//            System.out.println(key + "\t" + sameMass.get(key));
//          }
//          System.out.println("--------------------------");
//        }
//      }
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  // this is currently not required: rewrite if needed
//  private void translateTAGListToMzMineFormat()
//  {
//    try {
//      InputStream myxls = new FileInputStream("F:\\lipidomics\\Massenliste-TAG_mit IS.xls");
//      HSSFWorkbook workbook     = new HSSFWorkbook(myxls);
//      File fileForMzMine = new File("F:\\lipidomics\\20100210\\cdf\\masses.csv");
//      BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(fileForMzMine));
//      Hashtable<String,Vector<CgParameterSet>> sheetParams = new Hashtable<String,Vector<CgParameterSet>>();
//      long startCalcTime = System.currentTimeMillis();
//      for (int sheetNumber = 0; sheetNumber != workbook.getNumberOfSheets(); sheetNumber++) {
//        Sheet sheet = workbook.getSheet(sheetNumber);
//        int sideChainColumn = -1;
//        int doubleBondColumn = -1;
//        int carbonColumn = -1;
//        int hydrogenColumn = -1;
//        int deuteriumColumn = -1;
//        int oxygenColumn = -1;
//        int phosphateColumn = -1;
//        int nitrogenColumn = -1;
//        int massOfInterestColumn = -1;
//        int retTimeColumn = -1;
//        boolean foundColumns = false;
//        Vector<CgParameterSet> params = new Vector<CgParameterSet>();
//        float fixedStartTime = 0;
//        float fixedEndTime = Float.MAX_VALUE;
//        int idNumber = 1;
//        for (int rowCount = 0; rowCount != sheet.getRows(); rowCount++) {
//          String sideChain = "";
//          int doubleBonds = -1;
//          int carbon = -1;
//          int hydrogen = -1;
//          int deuterium = -1;
//          int oxygen = -1;
//          int phosphate = -1;
//          int nitrogen = -1;
//          float massOfInterest = -1;
//          float retTime = -1;
//          for (int i = 0; i != sheet.getColumns(); i++) {
//            String contents = sheet.getCell(i, rowCount).getContents();
//            if (contents != null)
//              contents = contents.trim();
//            if (!foundColumns) {
//              if (contents.equalsIgnoreCase("Seitenkette")) {
//                sideChainColumn = i;
//              }
//              if (contents.equalsIgnoreCase("dbs")
//                  || contents.equalsIgnoreCase("dbs_TAG")) {
//                doubleBondColumn = i;
//              }
//              if (contents.equalsIgnoreCase("C")) {
//                carbonColumn = i;
//              }
//              if (contents.equalsIgnoreCase("H")) {
//                hydrogenColumn = i;
//              }
//              if (contents.equalsIgnoreCase("D")) {
//                deuteriumColumn = i;
//              }
//              if (contents.equalsIgnoreCase("O")) {
//                oxygenColumn = i;
//              }
//              if (contents.equalsIgnoreCase("P")) {
//                phosphateColumn = i;
//              }
//              if (contents.equalsIgnoreCase("N")) {
//                nitrogenColumn = i;
//              }
//              if (contents.startsWith("mass") && contents.contains("(")
//                  && contents.contains(")")) {
//                massOfInterestColumn = i;
//              }
//              if (contents.equalsIgnoreCase("tR (min)")) {
//                retTimeColumn = i;
//              }
//              if (contents.startsWith("Start-RT:")) {
//                try {
//                  fixedStartTime = Float.parseFloat(contents.substring(
//                      "Start-RT:".length()).trim().replaceAll(",", "."));
//                }
//                catch (NumberFormatException nfx) {
//                  nfx.printStackTrace();
//                }
//                ;
//              }
//              if (contents.startsWith("Stop-RT:")) {
//                try {
//                  fixedEndTime = Float.parseFloat(contents.substring(
//                      "Stop-RT:".length()).trim().replaceAll(",", "."));
//                }
//                catch (NumberFormatException nfx) {
//                }
//                ;
//              }
//            } else {
//              if (i == sideChainColumn && contents != null
//                  & contents.length() > 0) {
//                sideChain = contents;
//              }
//              if (i == doubleBondColumn && contents != null
//                  && contents.length() > 0) {
//                doubleBonds = Integer.parseInt(contents);
//              }
//              if (i == carbonColumn && contents != null & contents.length() > 0) {
//                carbon = Integer.parseInt(contents);
//              }
//              if (i == hydrogenColumn && contents != null
//                  & contents.length() > 0) {
//                hydrogen = Integer.parseInt(contents);
//              }
//              if (i == deuteriumColumn && contents != null
//                  & contents.length() > 0) {
//                deuterium = Integer.parseInt(contents);
//              }
//              if (i == oxygenColumn && contents != null & contents.length() > 0) {
//                oxygen = Integer.parseInt(contents);
//              }
//              if (i == phosphateColumn && contents != null
//                  & contents.length() > 0) {
//                phosphate = Integer.parseInt(contents);
//              }
//              if (i == nitrogenColumn && contents != null
//                  & contents.length() > 0) {
//                nitrogen = Integer.parseInt(contents);
//              }
//              if (i == massOfInterestColumn && contents != null
//                  & contents.length() > 0) {
//                if (contents.contains(","))
//                  contents = contents.replaceAll(",", ".");
//                massOfInterest = Float.parseFloat(contents);
//              }
//              if (i == retTimeColumn && contents != null
//                  & contents.length() > 0) {
//                if (contents.contains(","))
//                  contents = contents.replaceAll(",", ".");
//                retTime = Float.parseFloat(contents);
//              }
//            }
//          }
//          if (sideChainColumn != -1 && carbonColumn != -1
//              && hydrogenColumn != -1 && oxygenColumn != -1
//              && phosphateColumn != -1 && nitrogenColumn != -1
//              && massOfInterestColumn != -1 && retTimeColumn != -1
//              && !foundColumns) {
//            foundColumns = true;
//          }
//          if (foundColumns && sideChain != null && sideChain.length() > 0
//              && massOfInterest > 0/* &&retTime>0 */) {
//            System.out.println("massOfInterest: " + massOfInterest);
//            String nameString = String.valueOf(sideChain);
//            if (doubleBonds > -1)
//              nameString += ":" + String.valueOf(doubleBonds);
//            nameString += "_";
//            String chemicalFormula = "";
//            if (carbon > 0) {
//              nameString += "C" + carbon;
//              if (chemicalFormula.length() > 0)
//                chemicalFormula += " ";
//              chemicalFormula += "C" + carbon;
//            }
//            if (hydrogen > 0) {
//              nameString += "H" + hydrogen;
//              if (chemicalFormula.length() > 0)
//                chemicalFormula += " ";
//              chemicalFormula += "H" + hydrogen;
//            }
//            if (deuterium > 0) {
//              nameString += "D" + deuterium;
//              if (chemicalFormula.length() > 0)
//                chemicalFormula += " ";
//              chemicalFormula += "D" + deuterium;
//            }
//            if (oxygen > 0) {
//              nameString += "O" + oxygen;
//              if (chemicalFormula.length() > 0)
//                chemicalFormula += " ";
//              chemicalFormula += "O" + oxygen;
//            }
//            if (phosphate > 0) {
//              nameString += "P" + phosphate;
//              if (chemicalFormula.length() > 0)
//                chemicalFormula += " ";
//              chemicalFormula += "P" + phosphate;
//            }
//            if (nitrogen > 0) {
//              nameString += "N" + nitrogen;
//              if (chemicalFormula.length() > 0)
//                chemicalFormula += " ";
//              chemicalFormula += "N" + nitrogen;
//            }
//            if (!nameString.substring(0,nameString.lastIndexOf("_")).equalsIgnoreCase("46:0")){
//
//              String toWrite = String.valueOf(idNumber)+"|"+String.valueOf(massOfInterest)+"|";
//              toWrite += "   |";
//              toWrite += "TG"+nameString.substring(0,nameString.lastIndexOf("_"))+"|"+nameString.substring(nameString.lastIndexOf("_")+1)+" (NH4)\n";
//              stream.write(toWrite.getBytes());
//              idNumber++;
//            }
//          }
//          
//        }
//        
//      }
//      stream.close();
//    }
//    catch (IOException e) {
//      // TODO Auto-generated catch block
//      e.printStackTrace();
//    }
//
//  }
  
  private void shortenMSList(){
    try{
      LineNumberReader reader = new LineNumberReader(new FileReader("F:\\lipidomics\\MSMSexamples\\PC-34-1-zu_20101125_MCT-KO-19_cleaned.txt"));
      String line = null;
      String output = "";
      Vector<String> mz = new Vector<String>();
      Vector<Double> intens = new Vector<Double>();
      double highest = 0;
      while ((line = reader.readLine())!=null){
        StringTokenizer tokenizer = new StringTokenizer(line,"\t");
        if (tokenizer.countTokens()==2){
          String mass = tokenizer.nextToken();
          String intString = tokenizer.nextToken();
          try{
            double intensity = Double.parseDouble(intString);
            if (intensity>500d){
              if (intensity>highest)
                highest = intensity;
              mz.add(mass);
              intens.add(intensity);
              //output += mass+"\t"+intString+"\n";
            }
          } catch (NumberFormatException nfx){
            
          }
        }
      }
      int count = 0;
      for (String mass : mz){
        output += mass+" "+String.valueOf(intens.get(count)/highest)+"\n";
        count++;
      }
      //System.out.println("("+output);
      reader.close();
      BufferedOutputStream stream  = new BufferedOutputStream(new FileOutputStream("F:\\lipidomics\\MSMSexamples\\PC-34-1-zu_20101125_MCT-KO-19_short.txt"));
      stream.write(output.getBytes());
      stream.close();
    }catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void calculateTheoreticalMass(){
    ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
    try {
      parser.parse();
      System.out.println(parser.calculateTheoreticalMass("C5 H11 O5 P1 N1", false));
      //System.out.println(parser.calculateTheoreticalMass("C44 H83 O8 P1 N1", false));
    }
    catch (SpectrummillParserException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  private void calculateByKnownMassesAndCAtoms(){
    try {
      ElementConfigParser aaParser = new ElementConfigParser("./elementNoCConfig.xml");
      aaParser.parse();
      ElementConfigParser probabilityParser = new ElementConfigParser("./elementNoIsotopesConfig.xml");
      probabilityParser.parse();
      InputStream myxls = new FileInputStream("F:\\TullnFH\\20111004\\helpDocuments\\found.xls");
      HSSFWorkbook workbook     = new HSSFWorkbook(myxls);
      Vector<KnowCAtomVO> vos = new Vector<KnowCAtomVO>();
      Vector<String> toSimilar = new  Vector<String>();
      for (int sheetNumber = 0; sheetNumber!=workbook.getNumberOfSheets(); sheetNumber++){
        HSSFSheet sheet = workbook.getSheetAt(sheetNumber);
        int numColumn = -1;
        int cColumn = -1;
        int massColumn = -1;
        int rtColumn = -1;
        String id = null;
        int cCount = -1;
        double mz = -1;
        double rt = -1;
        String previousId = null;
        double previousMass = 0;
        double previousRt = -1;
        int previousCAtoms = 0;
        boolean isPreviousActive = false;
        for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
          if (id!=null&&cCount>0&&mz>0){
            double rtDiff = rt-previousRt;
            if (rtDiff<0)rtDiff *= -1d;
            if ((mz-previousMass)>0.001f || previousCAtoms!=cCount || rtDiff>2){
//              if (isPreviousActive)vos.add(new KnowCAtomVO(previousId,previousCAtoms,previousMass,previousRt));
//              else 
              if (id.equals("29"))System.out.println(id+","+cCount+","+mz+","+rt);
              vos.add(new KnowCAtomVO(id,cCount,mz,rt));
              previousId = id;
              previousMass = mz;
              previousRt = rt;
              previousCAtoms = cCount;
//              isPreviousActive = false;
              previousId = id;
            }else{
              toSimilar.add(id);
//              isPreviousActive = true;
            }
              
          }
          id = null;
          cCount = -1;
          mz = -1;
          rt = -1;
          HSSFRow row = sheet.getRow(rowCount);
          for (int i=0; row!=null && i!=(row.getLastCellNum()+1);i++){
            HSSFCell cell = row.getCell(i);
            String contents = "";
            if (rowCount==0){
              if (cell==null)continue;
              contents = cell.getStringCellValue().trim();
              if (contents.equalsIgnoreCase("Num")) numColumn = i;
              if (contents.equalsIgnoreCase("CCount")) cColumn = i;
              if (contents.equalsIgnoreCase("mz")) massColumn = i;
              if (contents.equalsIgnoreCase("RT_min")) rtColumn = i;
//              System.out.println("numColumn: "+numColumn);
            }else{
              if (i==numColumn){
                id = String.valueOf((int)cell.getNumericCellValue());
              }else if (i==cColumn){
                cCount = (int)cell.getNumericCellValue();
              }else if (i == massColumn){
                mz = cell.getNumericCellValue();
              }else if (i == rtColumn){
                rt = cell.getNumericCellValue();
              }
            }
          }
        }
        if (id!=null&&cCount>0&&mz>0){
          double rtDiff = rt-previousRt;
          if (rtDiff<0)rtDiff *= -1d;
          if ((mz-previousMass)>0.001f || previousCAtoms!=cCount || rtDiff>2)
            vos.add(new KnowCAtomVO(id,cCount,mz,rt)); 
        }
      }
//      System.out.println("Vos: "+vos.size());
//      for (KnowCAtomVO vo:vos){
//        System.out.println(vo.getId());
//      }
      Hashtable<String,ChemicalFormulaVO> medianFormulas = new Hashtable<String,ChemicalFormulaVO>();
      Hashtable<String,String> usedElements = new Hashtable<String,String>();
      Vector<String> noFormulaFound = new Vector<String>();
      for (int i=0; i!=vos.size();i++){
        KnowCAtomVO vo = vos.get(i);
        System.out.println("Current one: "+vo.getId());
        double remainingMass = vo.getMz()-(double)(12*vo.getCCount());
        Vector<ChemicalFormulaVO> formulas = aaParser.calculatePossibleElementCompositions((float)remainingMass, 0.01f);
        if (formulas.size()>0){
//          System.out.println(i);
          Hashtable<String,Float> ratios = new Hashtable<String,Float>();
          Hashtable<String,ChemicalFormulaVO> formulaVOs = new Hashtable<String,ChemicalFormulaVO>();
          
          for (ChemicalFormulaVO formula : formulas){
            Vector<Double> distri = StaticUtils.calculateChemicalFormulaIntensityDistribution(probabilityParser, formula.getChemicalFormula(), 2, false);
            double zeroIsotope =  distri.get(0);
            double oneIsotope = 0;
            if (distri.size()>1)
              oneIsotope = distri.get(1);
            float ratio = (float)(oneIsotope/zeroIsotope);
//            System.out.println(formula.getChemicalFormula()+";"+ratio);
            ratios.put(formula.getChemicalFormula(),ratio);
            formulaVOs.put(formula.getChemicalFormula(), formula);
          }
          float median = Calculator.median(new Vector<Float>(ratios.values()));
          ChemicalFormulaVO medianFormula = null;
          float difference = Float.MAX_VALUE;
          for (String key: ratios.keySet()){
            float ratio = median/ratios.get(key);
            if (ratio<1) ratio = 1/ratio;
            float thisDiff = ratio-1;
            if (thisDiff<difference){
              difference=thisDiff;
              medianFormula = formulaVOs.get(key);
            }
          }
          medianFormulas.put(vo.getId(), medianFormula);
          for (SmChemicalElementVO elVO : medianFormula.getElements()){
            usedElements.put(elVO.getChemicalSymbol(), elVO.getChemicalSymbol());
          }
          
//        System.out.println(vo.getId()+";"+remainingMass+";"+formulas.size());
//        for (ChemicalFormulaVO formVO : formulas){
//          System.out.println(formVO.getChemicalFormula());
//        }
        }else{
          noFormulaFound.add(vo.getId());
        }
      }
      Vector<String> usedEles = new Vector<String>();
      usedEles.add("C");
      usedEles.add("Cc");
      for (String element : usedElements.keySet())usedEles.add(element);
      System.out.println("Finished calculating! Now storing");
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream("F:\\TullnFH\\20111004\\quantMet.xls"));
      BufferedOutputStream problematic = new BufferedOutputStream(new FileOutputStream("F:\\TullnFH\\20111004\\notPossible.txt"));
      HSSFWorkbook resultWorkbook = new HSSFWorkbook();
      HSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
      HSSFSheet resultSheet = resultWorkbook.createSheet("Metabolites");
      int currentRow = 3;
      HSSFRow row = resultSheet.createRow(currentRow);
      HSSFCell label = row.createCell(0,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("Name");
      label.setCellStyle(headerStyle);
      for (int i=0;i!=usedEles.size();i++){
        label = row.createCell((i+1),HSSFCell.CELL_TYPE_STRING);
        label.setCellValue(usedEles.get(i));
        label.setCellStyle(headerStyle);        
      }
      label = row.createCell(usedEles.size(),HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("mass(form[+H] name[])");
      label.setCellStyle(headerStyle);
      label = row.createCell(usedEles.size()+1,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("tR (min)");
      label.setCellStyle(headerStyle);      
      HSSFRow c13Row = null;
      for (int i=0; i!=vos.size();i++){
        KnowCAtomVO vo = vos.get(i);
        if (medianFormulas.containsKey(vo.getId())){
          ChemicalFormulaVO formVO = medianFormulas.get(vo.getId());
          currentRow++;
          row = resultSheet.createRow(currentRow);
          label = row.createCell(0,HSSFCell.CELL_TYPE_STRING);
          label.setCellValue(vo.getId());
          currentRow++;
          c13Row = resultSheet.createRow(currentRow);
          label = c13Row.createCell(0,HSSFCell.CELL_TYPE_STRING);
          label.setCellValue(vo.getId()+"-C13");
          for (int j=0;j!=usedEles.size();j++){
            int amount = 0;
            if (usedEles.get(j).equalsIgnoreCase("C")||usedEles.get(j).equalsIgnoreCase("Cc")){
              amount = vo.getCCount();
            }else{
              SmChemicalElementVO elem = new  SmChemicalElementVO(usedEles.get(j),new Vector<SmIsotopeVO>());
              amount = formVO.getNumberOfItemsOfOneElement(elem);
            }
            label = row.createCell((j+1),HSSFCell.CELL_TYPE_NUMERIC);
            if (!usedEles.get(j).equals("Cc")){
              label.setCellValue(amount);
            }else{
              label.setCellValue(0);
            }
            label = c13Row.createCell((j+1),HSSFCell.CELL_TYPE_NUMERIC);
            if (!usedEles.get(j).equals("C")){
              label.setCellValue(amount);
            }else{
              label.setCellValue(0);
            }              
          }
          label = row.createCell(usedEles.size(),HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(vo.getMz());
          label = row.createCell(usedEles.size()+1,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(vo.getRt());
          label = c13Row.createCell(usedEles.size(),HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(vo.getMz()+(1.003354826f*vo.getCCount()));
          label = c13Row.createCell(usedEles.size()+1,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(vo.getRt());

        }
      }
   
      resultWorkbook.write(out);
      out.close();
      String problemString = "";
      problemString+="To similar: \n";
      for (String similar: toSimilar){
        problemString+=similar+"\n";
      }
      problemString+="No formula found: \n";
      for (String notFormula : noFormulaFound){
        problemString+=notFormula+"\n";
      }
      problematic.write(problemString.getBytes());
      problematic.close();
    }
    catch (SpectrummillParserException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }catch (IOException e) {
      e.printStackTrace();
//    } catch (WriteException e) {
//      // TODO Auto-generated catch block
//      e.printStackTrace();
    }
    
  }

  private class KnowCAtomVO{
    
    private String id_;
    private int cCount_;
    private double mz_;
    private double rt_;
    
    public KnowCAtomVO(String id, int cCount, double mz, double rt)
    {
      id_ = id;
      cCount_ = cCount;
      mz_ = mz;
      rt_ = rt;
    }

    public String getId()
    {
      return id_;
    }

    public int getCCount()
    {
      return cCount_;
    }

    public double getMz()
    {
      return mz_;
    }
    
    public double getRt()
    {
      return rt_;
    }
    
  }
  
  private static HSSFCellStyle getHeaderStyle(HSSFWorkbook wb){
    HSSFCellStyle arial12style = wb.createCellStyle();
    HSSFFont arial12font = wb.createFont();
    arial12font.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(HSSFCellStyle.ALIGN_CENTER);
    return arial12style;
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
  
  private void calculateIntensityDistribution(){
	
    ElementConfigParser aaParser = Settings.getElementParser();
    try {
//      Vector<Double> mustMatchProbabs = aaParser.calculateChemicalFormulaIntensityDistribution("H22 P1 O1 Cc15", 2, false);
//      Vector<Double> mustMatchProbabs = aaParser.calculateChemicalFormulaIntensityDistribution("Cc33 H59 S2 F6 P2 O2", 2, false,true);
      //Vector<Double> mustMatchProbabs = aaParser.calculateChemicalFormulaIntensityDistribution("Cc10 H18 P1 O2", 5, false).get(0);
      //Vector<Double> mustMatchProbabs = aaParser.calculateChemicalFormulaIntensityDistribution("C6 H12 O6", 5, false).get(0);
//      Vector<Double> mustMatchProbabs = StaticUtils.calculateChemicalFormulaIntensityDistribution(aaParser, "C100 O26 H95", 5, false);
      //Vector<Double> mustMatchProbabs = StaticUtils.calculateChemicalFormulaIntensityDistribution(aaParser, "Cc6 H11 O6", 5, false);
      Vector<Double> mustMatchProbabs = StaticUtils.calculateChemicalFormulaIntensityDistribution(aaParser, "C6 H12 O6", 5, false);
      
      for (Double prob : mustMatchProbabs){
        System.out.println(prob);
      }
    }
    catch (SpectrummillParserException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

  }
  
  private void justMzValues(){
    try{
      LineNumberReader reader = new LineNumberReader(new FileReader("E:\\lipidomics\\MSMSexamples\\PE_36-4-zu_20101125_MCT-KO-19_cleaned.txt"));
      BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream("E:\\lipidomics\\MSMSexamples\\PE_36-4-zu_20101125_MCT-KO-19_cleanedJustMz.txt"));
      String line;
      Vector<String> mz = new Vector<String>();
      Vector<String> intens = new Vector<String>();
      double highestInt = 0;
      while ((line = reader.readLine()) != null) {
        StringTokenizer tokenizer = new StringTokenizer(line,"\t ");
        String newFormat = "";
        newFormat = tokenizer.nextToken();
        newFormat = newFormat.substring(0,newFormat.indexOf(".")+4);
        mz.add(newFormat);
//        newFormat += " ";
        newFormat = tokenizer.nextToken();
        intens.add(newFormat);
        double intensity = Double.parseDouble(newFormat);
        if (intensity>highestInt) highestInt = intensity;
//        newFormat = newFormat.substring(0,newFormat.indexOf(".")+4);
//        newFormat += "\n";
//        stream.write(newFormat.getBytes());
      }
      for (int i=0; i!=mz.size();i++){
        double relInt = (Double.parseDouble(intens.get(i))*100d)/highestInt;
        String relIntString = String.valueOf(Calculator.roundDBL(relInt, 3));
        if (relIntString.indexOf(".")==-1&&relIntString.length()>(relIntString.indexOf(".")+4))
          relIntString = relIntString.substring(0,relIntString.indexOf(".")+4);
        String newLine = mz.get(i)+" "+relIntString+"\n";
        stream.write(newLine.getBytes());
      }
      reader.close();
      stream.close();
    }catch(Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void ttest(){
    int sampleSize1=14;
    int sampleSize2=21;
    double mean1 = 12.40357143d;
    double stdev1 = 1.415019575d;
    double mean2 = 70.79333333d;
    double stdev2 = 5.915576824d;  
////    IndependentTwoSamplesTTest ttest = new IndependentTwoSamplesTTest(sampleSize1,sampleSize2, mean1, mean2, stdev1, stdev2);
////    double[] results = ttest.getResults();
////    System.out.println("p-Value: "+results[0]+"; df-double: "+results[1]+"; df: "+results[2]);
  }
  
  private void ttestFromExcel(){
    try {
//      InputStream myxls = new FileInputStream("L:\\HFD-Paper\\heatmaps\\stdevs\\PI_n3.xls");
      InputStream myxls = new FileInputStream("E:\\lipidomics\\20120823\\ATGL-p-values.xls");
      HSSFWorkbook workbook     = new HSSFWorkbook(myxls);
      HSSFSheet sheet = workbook.getSheetAt(0);
      HSSFRow row = sheet.getRow(1);
      Vector<Integer> columnsOfInterest = new Vector<Integer>();
      Hashtable<Integer,String> names = new Hashtable<Integer,String>();
      for (int i=0;  row!=null && i!=(row.getLastCellNum()+1);i++){
        HSSFCell cell = row.getCell(i);
        String contents = "";
        int cellType = -1;
        if (cell!=null) cellType = cell.getCellType();
        if (cellType==HSSFCell.CELL_TYPE_STRING){
          contents = cell.getStringCellValue();
        }else if (cellType==HSSFCell.CELL_TYPE_NUMERIC || cellType==HSSFCell.CELL_TYPE_FORMULA){
         double numeric = cell.getNumericCellValue();
         contents = String.valueOf(numeric);
        }
//        if (contents.length()>2 && contents.toCharArray()[2]==':'){
        if (contents.length()==2){
          columnsOfInterest.add(i);
          names.put(i, contents);
        }
      }
      Hashtable<Integer,Double> fedMeans = readValuesOfTTestRow(sheet.getRow(2),names);
      Hashtable<Integer,Double> fedStdevs = readValuesOfTTestRow(sheet.getRow(3),names);
      Hashtable<Integer,Double> fasMeans = readValuesOfTTestRow(sheet.getRow(4),names);
      Hashtable<Integer,Double> fasStdevs = readValuesOfTTestRow(sheet.getRow(5),names);
      Hashtable<Integer,Double> hfdMeans = readValuesOfTTestRow(sheet.getRow(6),names);
      Hashtable<Integer,Double> hfdStdevs = readValuesOfTTestRow(sheet.getRow(7),names);
      Hashtable<Integer,Double> koFasMeans = readValuesOfTTestRow(sheet.getRow(8),names);
      Hashtable<Integer,Double> koFasStdevs = readValuesOfTTestRow(sheet.getRow(9),names);

      for (Integer i:columnsOfInterest){
        String name = names.get(i);
        System.out.println(name+":");
        this.printResults("KO-FED/WT-FED", fedMeans.get(i), fasMeans.get(i), fedStdevs.get(i), fasStdevs.get(i));
        this.printResults("WT-FAS/WT-FED", fedMeans.get(i), hfdMeans.get(i), fedStdevs.get(i), hfdStdevs.get(i));
        this.printResults("KO-FAS/WT-FED", fedMeans.get(i), koFasMeans.get(i), fedStdevs.get(i), koFasStdevs.get(i));
        this.printResults("WT-FAS/KO-FED", fasMeans.get(i), hfdMeans.get(i), fasStdevs.get(i), hfdStdevs.get(i));
        this.printResults("KO-FAS/KO-FED", fasMeans.get(i), koFasMeans.get(i), fasStdevs.get(i), koFasStdevs.get(i));
        this.printResults("WT-FAS/KO-FAS", hfdMeans.get(i), koFasMeans.get(i), hfdStdevs.get(i), koFasStdevs.get(i));

        
      }
      myxls.close();
    }catch(Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void printResults(String startString, double mean1, double mean2, double stdev1, double stdev2){
    int sampleSize = 3;
//    double[] results = (new IndependentTwoSamplesTTest(sampleSize,sampleSize, mean1, mean2, stdev1, stdev2)).getResults();
//    System.out.println(startString+": "+results[0]+" ; "+results[1]+" ; "+results[2]);
  }
  
  private Hashtable<Integer,Double> readValuesOfTTestRow(HSSFRow row, Hashtable<Integer,String> names){
    Hashtable<Integer,Double> vals = new Hashtable<Integer,Double>();
    for (int i=0;  row!=null && i!=(row.getLastCellNum()+1);i++){
      if (names.containsKey(i)){
        HSSFCell cell = row.getCell(i);
        Double numeric = null;
        int cellType = -1;
        if (cell!=null) cellType = cell.getCellType();
        if (cellType==HSSFCell.CELL_TYPE_STRING){
          String contents = cell.getStringCellValue();
          try{ numeric = new Double(contents);}catch(NumberFormatException nfx){};
        }else if (cellType==HSSFCell.CELL_TYPE_NUMERIC || cellType==HSSFCell.CELL_TYPE_FORMULA){
         numeric = cell.getNumericCellValue();
        }
        vals.put(i, numeric);
      }
    }
    return vals;
  }
  
  private void msFileReaderAccess(){
    System.loadLibrary("XRawfile2_x64");
    System.loadLibrary("fregistry_x64");
    System.loadLibrary("Fileio_x64");
    
    //XRawfileCtrl ctrl = XRawfileCtrl.INSTANCE;
//    NativeLibrary.addSearchPath("XRawfile2_x64", "./XRawfile2_x64.dll");
//    NativeLibrary.addSearchPath("fregistry_x64", "./fregistry_x64.dll");
//    NativeLibrary.addSearchPath("Fileio_x64", "./Fileio_x64.dll");
//    IXRawfile ctrl = at.tugraz.genome.thermo.ClassFactory.createMSFileReader_XRawfile();    
//    System.out.println("I have an instance");
//    
//    ctrl.open("E:\\lipidomics\\20100210\\20100126_TAG-34.RAW");
//    ctrl.setCurrentController(0,1);
//    System.out.println("The file is open");
//    Holder<Integer> stringHolder = new Holder<Integer>();
//    ctrl.getLastSpectrumNumber(stringHolder);
//    System.out.println(stringHolder.value);
//    ctrl.close();
//    System.out.println("The file is closed");
  }
  

  private void paintSpectrumFromXls(){
    try{
      InputStream myxls = new FileInputStream("L:\\FWFAntrag\\fragmentierungsBeispiele\\TG_878.xls");
      HSSFWorkbook workbook     = new HSSFWorkbook(myxls);
      Vector<Double> masses = new Vector<Double>();
      Vector<Double> intensities = new Vector<Double>();
      double highestIntensity=0d;
      for (int sheetNumber = 0; sheetNumber!=workbook.getNumberOfSheets(); sheetNumber++){
        HSSFSheet sheet = workbook.getSheetAt(sheetNumber); 
        int massColumn = -1;
        int intColumn = -1;
        for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
          HSSFRow row = sheet.getRow(rowCount);
          double mass = -1d;
          double intensity = -1d;
          for (int i=0; row!=null && i!=(row.getLastCellNum()+1);i++){
            HSSFCell cell = row.getCell(i);
            String contents = "";
            Double numeric = null;
            int cellType = -1;
            if (cell!=null) cellType = cell.getCellType();
            if (cellType==HSSFCell.CELL_TYPE_STRING){
              contents = cell.getStringCellValue();
              try{ numeric = new Double(contents);}catch(NumberFormatException nfx){};
            }else if (cellType==HSSFCell.CELL_TYPE_NUMERIC || cellType==HSSFCell.CELL_TYPE_FORMULA){
             numeric = cell.getNumericCellValue();
             contents = String.valueOf(numeric);
            }
            //String contents = sheet.getCell(i,rowCount).getContents();
            if (contents!=null)
              contents = contents.trim();
            if (massColumn>=0 && intColumn>=0 && numeric!=null){
              if (i==massColumn) mass = numeric;
              if (i==intColumn) intensity = numeric;
            }else{
              if (contents.equalsIgnoreCase("Mass")) massColumn=i;
              if (contents.equalsIgnoreCase("Intensity")) intColumn=i;
            }
          }
          if (mass>-0.5d && intensity>-0.5d){
            masses.add(mass);
            intensities.add(intensity);
            if (intensity>highestIntensity) highestIntensity=intensity;
          }
        }
      }
      myxls.close();
//      XYSeries series = new XYSeries("spectrum");
//      for (int i=0; i!=masses.size(); i++){
//        Double intensity = (intensities.get(i)*100d)/highestIntensity;
//        if (intensity>0.0d) series.add(masses.get(i), intensity);
//      }
//      XYSeriesCollection coll = new XYSeriesCollection(series);
//      JFreeChart chart = ChartFactory.createXYBarChart("Scan 878", "m/z",
//          false, "Relative Intensity [%]", coll, PlotOrientation.VERTICAL,
//          true, true, false);
//      XYPlot plot = chart.getXYPlot(); 
//      XYBarRenderer renderer = new XYBarRenderer();
//      renderer.setSeriesPaint(0,Color.BLACK);
//      plot.setRenderer(renderer);
//
//      BufferedOutputStream stream = new BufferedOutputStream(
//          new FileOutputStream("L:\\FWFAntrag\\fragmentierungsBeispiele\\TG_878.png"));
//      ImageIO.write(chart.createBufferedImage(1000, 700), "PNG",
//          stream);
//      stream.close();

    } catch (Exception ex){
      
    }
  }

  private void plus2IsotopicRelationXls(){
    try{
//      Vector contents = QuantificationThread.parseQuantExcelFile("E:\\lipidomics\\Niki\\20140515_BMP_pos.xls", 2f, 3f, 3, 1, true, 0.0001f, 0f, -1f, -1f);
//      LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)contents.get(0);
//      Hashtable<String,Vector<String>> analyteSequence = (Hashtable<String,Vector<String>>)contents.get(1);
//      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)contents.get(2);
//      
//      
//      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream("E:\\lipidomics\\20121108\\isoDistribution.xls"));   
//      HSSFWorkbook resultWorkbook = new HSSFWorkbook();
//      HSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
//      for (String lclass : analyteSequence.keySet()){
//        Vector<String> analyteNames = analyteSequence.get(lclass);        
//        HSSFSheet resultSheet = resultWorkbook.createSheet(lclass);
//        
//        HSSFCell number;
//        HSSFRow row = resultSheet.createRow(0);
//        HSSFCell label = row.createCell(0,HSSFCell.CELL_TYPE_STRING);
//        label.setCellValue("Name");
//        label.setCellStyle(headerStyle);
//        label = row.createCell(1,HSSFCell.CELL_TYPE_STRING);
//        label.setCellValue("Iso +0");
//        label.setCellStyle(headerStyle);    
//        label = row.createCell(2,HSSFCell.CELL_TYPE_STRING);
//        label.setCellValue("Iso +1");
//        label.setCellStyle(headerStyle);
//        label = row.createCell(3,HSSFCell.CELL_TYPE_STRING);
//        label.setCellValue("Iso +2");
//        label.setCellStyle(headerStyle);
//        label = row.createCell(4,HSSFCell.CELL_TYPE_STRING);
//        label.setCellValue("Iso +3");
//        label.setCellStyle(headerStyle);
//        for (int i=0; i!=analyteNames.size();i++){
//          String analyte = analyteNames.get(i);
//          Hashtable<String,QuantVO> modHash = quantObjects.get(lclass).get(analyte);
//          int amountOfModifications = modHash.size();
//          String name = new String(analyte);
//          Set<String> mods = modHash.keySet();
//          int count=0;
//          for (String mod : mods ){
//            QuantVO quant = modHash.get(mod);
//            Vector<Double> probabs = quant.getProbabs();
//            if (amountOfModifications>1) name+="_"+mod;
//            HSSFRow dataRow = resultSheet.createRow(i*amountOfModifications+1+count);
//            label = dataRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
//            label.setCellValue(name);
//            number = dataRow.createCell(1,HSSFCell.CELL_TYPE_NUMERIC);
//            number.setCellValue(probabs.get(0));
//            number = dataRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
//            number.setCellValue(probabs.get(1));
//            number = dataRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
//            number.setCellValue(probabs.get(2));
//            number = dataRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
//            number.setCellValue(probabs.get(3));
//
//            count++;  
//          }
//        }
//      }
//      resultWorkbook.write(out);
//      out.close();
//
    }catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void readQuantExcel(){
    try{
////      Vector excelContent = QuantificationThread.parseQuantExcelFile("D:\\Kim\\SL_massLists\\SL_negative.xlsx", -1f, -1f, 2, 1, true, 0f, 0f, -1f, -1f);
      Vector excelContent = null;
      System.out.println("Hallo");
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)excelContent.get(3);
      LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)excelContent.get(0);
      Hashtable<String,Vector<String>> analyteSequence = (Hashtable<String,Vector<String>>)excelContent.get(1);

      for (String className :classSequence.keySet()){
        Hashtable<String,Hashtable<String,QuantVO>> classQuant = quantObjects.get(className);
        Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultsClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
        System.out.println("----------------- "+className+" -------------------");
//        if (className.equalsIgnoreCase("oxPC_38_4")||className.equalsIgnoreCase("oxPC_40_6")){
          for (String analyteName : analyteSequence.get(className)){
            Hashtable<String,QuantVO> analyteQuant = classQuant.get(analyteName);
            System.out.println(analyteName);
            //System.out.println(analyteQuant.keySet());
//            if (analyteName.equalsIgnoreCase("SONPC")||analyteName.equalsIgnoreCase("7-OH-10O-4,8decendienoyl")){
//              System.out.println(analyteName);
//              QuantVO vo = analyteQuant.values().iterator().next();
//              System.out.println(vo.getAnalyteFormula());
//            }
          }  
//        }
      }  
    }catch(Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void translateMarleneFileToQuant(){
    try{
      ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
      parser.parse();
      String marleneFile = "E:\\Marlene\\Clematis_compounds_17122012.xlsx";
////      String quantFile = "E:\\Marlene\\Clematis_compounds.xls";
      InputStream myxls = new FileInputStream(marleneFile);
      XSSFWorkbook workbook     = new XSSFWorkbook(myxls);
      XSSFSheet sheet = workbook.getSheetAt(0);
      
////      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(quantFile));
      Hashtable<String,String> names = new Hashtable<String,String>();
      HSSFWorkbook resultWorkbook = new HSSFWorkbook();
      HSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
      HSSFSheet resultSheet = resultWorkbook.createSheet("Clemantis");
      int rowCountOut = 3;
      HSSFRow outRow = resultSheet.createRow(rowCountOut);
      HSSFCell label = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("Name");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("C");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("H");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("O");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("N");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("form[-H] name[-H]");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("tR (min)");
      label.setCellStyle(headerStyle);
      rowCountOut++;
      
      for (int rowCount=1;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
        XSSFRow row = sheet.getRow(rowCount);
        outRow = resultSheet.createRow(rowCountOut);
        
        XSSFCell cell = row.getCell(0);
        String name = cell.getStringCellValue();
        
        String formulaString = "";
        cell = row.getCell(3);
        String cString = cell.getStringCellValue();
        int amountC = 0;
        cell = row.getCell(4);
        Double value = getNumericCellValue(cell);
        if (value!=null) amountC = (int)Math.round(value);
        if (amountC > 0)
          formulaString += cString+String.valueOf(amountC);
        cell = row.getCell(5);
        String hString = cell.getStringCellValue();
        int amountH = 0;
        cell = row.getCell(6);
        value = getNumericCellValue(cell);
        if (value!=null) amountH = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+hString+String.valueOf(amountH);
        cell = row.getCell(7);
        String oString = cell.getStringCellValue();
        int amountO = 0;
        cell = row.getCell(8);
        value = getNumericCellValue(cell);
        if (value!=null) amountO = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+oString+String.valueOf(amountO);
        cell = row.getCell(9);
        String nString = cell.getStringCellValue();
        int amountN = 0;
        cell = row.getCell(10);
        value = getNumericCellValue(cell);
        if (value!=null) amountN = (int)Math.round(value);
        if (amountN > 0)
          formulaString += " "+nString+String.valueOf(amountN);
        
        if (formulaString.length()>0){
          HSSFCell cellOut = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
          cellOut.setCellValue(name);
          if (names.containsKey(name)) System.out.println("Duplicate: "+name);
          names.put(name, name);
          cellOut = outRow.createCell(1,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(amountC);
          cellOut = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(amountH);
          cellOut = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(amountO);
          cellOut = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(amountN);
          cellOut = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(parser.calculateTheoreticalMass(formulaString, false));
        }
        rowCountOut++;
      }
////      resultWorkbook.write(out);
////      out.close();
      myxls.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }
  }
    
  private Double getNumericCellValue(XSSFCell cell){
    String contents = "";
    Double numeric = null;   
    int cellType = -1;
    if (cell!=null) cellType = cell.getCellType();
    if (cellType==HSSFCell.CELL_TYPE_STRING){
      contents = cell.getStringCellValue();
      try{ numeric = new Double(contents);}catch(NumberFormatException nfx){};
    }else if (cellType==HSSFCell.CELL_TYPE_NUMERIC || cellType==HSSFCell.CELL_TYPE_FORMULA){
      numeric = cell.getNumericCellValue();
    }
    return numeric;
  }
  
  private Double getNumericCellValue(HSSFCell cell){
    String contents = "";
    Double numeric = null;   
    int cellType = -1;
    if (cell!=null) cellType = cell.getCellType();
    if (cellType==HSSFCell.CELL_TYPE_STRING){
      contents = cell.getStringCellValue();
      try{ numeric = new Double(contents);}catch(NumberFormatException nfx){};
    }else if (cellType==HSSFCell.CELL_TYPE_NUMERIC || cellType==HSSFCell.CELL_TYPE_FORMULA){
      numeric = cell.getNumericCellValue();
    }
    return numeric;
  }

  
  private void mergeExcelResults(){
    File dir = new File("E:\\Marlene\\AuswertungAndreas");
    File[] files = dir.listFiles();
    Hashtable<String,String> fileStarts = new Hashtable<String,String>();
    for (File file : files){
      String absPath = file.getAbsolutePath();
      String uniqueRawName = null;
      if (absPath.endsWith("_Clemantis_compounds.xlsx")) uniqueRawName = absPath.substring(0,absPath.length()-("_Clemantis_compounds.xlsx".length()));
      if (absPath.endsWith("_Clemantis_compounds_100-300.xlsx")) uniqueRawName = absPath.substring(0,absPath.length()-("_Clemantis_compounds_100-300.xlsx".length()));
      if (uniqueRawName!=null)
        fileStarts.put(uniqueRawName, uniqueRawName);
    }
    System.out.println(fileStarts.size());
    for (String fileStart : fileStarts.keySet()){ 
    String filePath1 = fileStart+"_Clemantis_compounds_100-300.xlsx";
    String filePath2 = fileStart+"_Clemantis_compounds_300-500.xlsx";
    String filePath3 = fileStart+"_Clemantis_compounds_500-700.xlsx";
    String filePath4 = fileStart+"_Clemantis_compounds_700-900.xlsx";
    String filePath5 = fileStart+"_Clemantis_compounds_900-1100.xlsx";
    String filePath6 = fileStart+"_Clemantis_compounds_1100-1300.xlsx";
    String filePath7 = fileStart+"_Clemantis_compounds_1300-4000.xlsx";
    Hashtable<String,Boolean> resultsShowModification1 = new Hashtable<String,Boolean>();
    Hashtable<String,Boolean> resultsShowModification2 = new Hashtable<String,Boolean>();
    Hashtable<String,Boolean> resultsShowModification3 = new Hashtable<String,Boolean>();
    Hashtable<String,Boolean> resultsShowModification4 = new Hashtable<String,Boolean>();
    Hashtable<String,Boolean> resultsShowModification5 = new Hashtable<String,Boolean>();
    Hashtable<String,Boolean> resultsShowModification6 = new Hashtable<String,Boolean>();
    Hashtable<String,Boolean> resultsShowModification7 = new Hashtable<String,Boolean>();
    try {
      QuantificationResult result1 = LDAResultReader.readResultFile(filePath1,  resultsShowModification1);
//      Hashtable<String,Vector<LipidParameterSet>> resultParams2 = new Hashtable<String,Vector<LipidParameterSet>>();
//      Hashtable<String,Integer> msLevels2 = new Hashtable<String,Integer>();
      QuantificationResult result2 = LDAResultReader.readResultFile(filePath2,  resultsShowModification2);
      QuantificationResult result3 = LDAResultReader.readResultFile(filePath3,  resultsShowModification3);
      QuantificationResult result4 = LDAResultReader.readResultFile(filePath4,  resultsShowModification4);
      QuantificationResult result5 = LDAResultReader.readResultFile(filePath5,  resultsShowModification5);
      QuantificationResult result6 = LDAResultReader.readResultFile(filePath6,  resultsShowModification6);
      QuantificationResult result7 = LDAResultReader.readResultFile(filePath7,  resultsShowModification7);
      Hashtable<String,Vector<String>> correctAnalyteSequence = (Hashtable<String,Vector<String>>)QuantificationThread.getCorrectAnalyteSequence("E:\\Marlene\\201403_untargeted\\Clemantis_compounds.xls",false).get(1);
      
      QuantificationResult mergedResults = mergeResults(result1,result2,correctAnalyteSequence);
      mergedResults = mergeResults(mergedResults,result3,correctAnalyteSequence);
      mergedResults = mergeResults(mergedResults,result4,correctAnalyteSequence);
      mergedResults = mergeResults(mergedResults,result5,correctAnalyteSequence);
      mergedResults = mergeResults(mergedResults,result6,correctAnalyteSequence);
      mergedResults = mergeResults(mergedResults,result7,correctAnalyteSequence);
      
//      Hashtable<String,Vector<LipidParameterSet>> mergedParams = new Hashtable<String,Vector<LipidParameterSet>>();
//      Map<String,Integer> mergedLevels = null;
//      
//      for (String lClass : correctAnalyteSequence.keySet()){
//        Vector<String> anals = correctAnalyteSequence.get(lClass);
//        Vector<LipidParameterSet> anals1 = result1.getIdentifications().get(lClass);
//        Vector<LipidParameterSet> anals2 = result2.getIdentifications().get(lClass);
//        Vector<LipidParameterSet> analsMerged = new Vector<LipidParameterSet>();
//        for (String analName : anals){
//          Vector<LipidParameterSet> sameAnalsDiffRt = new Vector<LipidParameterSet>();
//          for (LipidParameterSet set : anals1){
//            if (analName.equalsIgnoreCase(set.getName())){
//              int position = -1;
//              for (int i=(sameAnalsDiffRt.size()-1);i!=-1;i--){
//                if ((Double.parseDouble(set.getRt())+0.001d)<Double.parseDouble(sameAnalsDiffRt.get(i).getRt())) position = i;
//              }
//              if (position == -1) sameAnalsDiffRt.add(set);
//              else sameAnalsDiffRt.add(position,set);
//            }
//          }
//          for (LipidParameterSet set : anals2){
//            if (analName.equalsIgnoreCase(set.getName())){
//              int position = -1;
//              for (int i=(sameAnalsDiffRt.size()-1);i!=-1;i--){
//                if ((Double.parseDouble(set.getRt())+0.001d)<Double.parseDouble(sameAnalsDiffRt.get(i).getRt())) position = i;
//              }
//              if (position == -1) sameAnalsDiffRt.add(set);
//              else sameAnalsDiffRt.add(position,set);
//            }
//          }
//          for (LipidParameterSet set : sameAnalsDiffRt) analsMerged.add(set);
//        }
//        mergedParams.put(lClass, analsMerged);
//      }
      //String resultFile = fileStart+"_Clemantis_compounds.xlsx";
//      QuantificationResult mergedResults = new QuantificationResult(mergedParams,result1.getConstants(),mergedLevels);
// To prevent overwrite - commented      
      int directoryIndex = fileStart.lastIndexOf("/");
      if (fileStart.lastIndexOf("\\")>directoryIndex) directoryIndex = fileStart.lastIndexOf("\\");
      String resultFile = fileStart.substring(0,directoryIndex)+File.separator+"merged"+fileStart.substring(directoryIndex)+"_Clemantis_compounds.xlsx";
      QuantificationThread.writeResultsToExcel(resultFile,mergedResults);
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    }
  }
  
  private QuantificationResult mergeResults (QuantificationResult result1, QuantificationResult result2, Hashtable<String,Vector<String>> correctAnalyteSequence){
    Hashtable<String,Vector<LipidParameterSet>> mergedParams = new Hashtable<String,Vector<LipidParameterSet>>();
    Map<String,Integer> mergedLevels = null;
    
    for (String lClass : correctAnalyteSequence.keySet()){
      Vector<String> anals = correctAnalyteSequence.get(lClass);
      Vector<LipidParameterSet> anals1 = result1.getIdentifications().get(lClass);
      Vector<LipidParameterSet> anals2 = result2.getIdentifications().get(lClass);
      Vector<LipidParameterSet> analsMerged = new Vector<LipidParameterSet>();
      for (String analName : anals){
        Vector<LipidParameterSet> sameAnalsDiffRt = new Vector<LipidParameterSet>();
        for (LipidParameterSet set : anals1){
          if (analName.equalsIgnoreCase(set.getName())){
            int position = -1;
            for (int i=(sameAnalsDiffRt.size()-1);i!=-1;i--){
              if ((Double.parseDouble(set.getRt())+0.001d)<Double.parseDouble(sameAnalsDiffRt.get(i).getRt())) position = i;
            }
            if (position == -1) sameAnalsDiffRt.add(set);
            else sameAnalsDiffRt.add(position,set);
          }
        }
        for (LipidParameterSet set : anals2){
          if (analName.equalsIgnoreCase(set.getName())){
            int position = -1;
            for (int i=(sameAnalsDiffRt.size()-1);i!=-1;i--){
              if ((Double.parseDouble(set.getRt())+0.001d)<Double.parseDouble(sameAnalsDiffRt.get(i).getRt())) position = i;
            }
            if (position == -1) sameAnalsDiffRt.add(set);
            else if (position>0 && sameAnalsDiffRt.get(position-1).getRt().equalsIgnoreCase(set.getRt())) sameAnalsDiffRt.add(position-1,set);
            else sameAnalsDiffRt.add(position,set);
          }
        }
        for (LipidParameterSet set : sameAnalsDiffRt) analsMerged.add(set);
      }
      mergedParams.put(lClass, analsMerged);
    }
    QuantificationResult mergedResults = new QuantificationResult(mergedParams,result1.getConstants(),mergedLevels,Settings.getFaHydroxyEncoding(),
        Settings.getLcbHydroxyEncoding());
    return mergedResults;
  }
  
  private void calcPotentialFormulasPerMass(){
    ElementConfigParser aaParser = new ElementConfigParser("./elementNoIsotopesConfig.xml");
    try {
      aaParser.parse();
      System.out.println("Hallo");
      Vector<ChemicalFormulaVO> formulas = aaParser.calculatePossibleElementCompositions(173.09900f, 0.005f);
      //173.0086
      for (ChemicalFormulaVO formula : formulas){
        boolean print = false;
        for (SmChemicalElementVO smVO : formula.getElements()){
          if (smVO.getChemicalSymbol().equalsIgnoreCase("C") && formula.getNumberOfItemsOfOneElement(smVO)>3){
            print = true;
          }    
        }
        if (print)
          System.out.println(formula.getChemicalFormula());
      }
    }
    catch (SpectrummillParserException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  private void translateSabineFileToQuant(){
    try{
      ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
      parser.parse();
      String marleneFile = "E:\\Sabine\\Lonicera_compounds_corrected.xlsx";
      String quantFile = "E:\\Sabine\\Lonicera_compounds_named.xls";
      InputStream myxls = new FileInputStream(marleneFile);
      XSSFWorkbook workbook     = new XSSFWorkbook(myxls);
      XSSFSheet sheet = workbook.getSheetAt(0);
      
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(quantFile));
      Hashtable<String,String> names = new Hashtable<String,String>();
      Hashtable<String,String> formula = new Hashtable<String,String>();
      HSSFWorkbook resultWorkbook = new HSSFWorkbook();
      HSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
      HSSFSheet resultSheet = resultWorkbook.createSheet("Lonicera");
      int rowCountOut = 3;
      HSSFRow outRow = resultSheet.createRow(rowCountOut);
      HSSFCell label = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("Name");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("C");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("H");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("O");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("N");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("Cl");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("form[-H] name[-H]");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(8,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("tR (min)");
      label.setCellStyle(headerStyle);
      rowCountOut++;
      //Juergen: attention - the last row number is fixed
      for (int rowCount=2;rowCount!=(104);rowCount++){
        XSSFRow row = sheet.getRow(rowCount);
        outRow = resultSheet.createRow(rowCountOut);
//        System.out.println(rowCount);
        XSSFCell cell = row.getCell(0);
        String name = cell.getStringCellValue();
        String formulaString = "";
        
        String cString = "C";
        int amountC = 0;
        cell = row.getCell(4);
        Double value = getNumericCellValue(cell);
        if (value!=null) amountC = (int)Math.round(value);
        if (amountC > 0)
          formulaString += cString+String.valueOf(amountC);
        String hString = "H";
        int amountH = 0;
        cell = row.getCell(5);
        value = getNumericCellValue(cell);
        if (value!=null) amountH = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+hString+String.valueOf(amountH);
        String oString = "O";
        int amountO = 0;
        cell = row.getCell(6);
        value = getNumericCellValue(cell);
        if (value!=null) amountO = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+oString+String.valueOf(amountO);
        String nString = "N";
        int amountN = 0;
        cell = row.getCell(7);
        value = getNumericCellValue(cell);
        if (value!=null) amountN = (int)Math.round(value);
        if (amountN > 0)
          formulaString += " "+nString+String.valueOf(amountN);
        String clString = "Cl";
        int amountCl = 0;
        cell = row.getCell(8);
        value = getNumericCellValue(cell);
        if (value!=null) amountCl = (int)Math.round(value);
        if (amountCl > 0)
          formulaString += " "+clString+String.valueOf(amountCl);

        
        if (formulaString.length()>0){
          HSSFCell cellOut = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
          cellOut.setCellValue(name);
          if (names.containsKey(name)) System.out.println("Duplicate: "+name);
          names.put(name, name);
          if (formula.containsKey(formulaString)) System.out.println("Duplicate: "+name+" ; "+formulaString);
          formula.put(formulaString, formulaString);

          cellOut = outRow.createCell(1,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(amountC);
          cellOut = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(amountH);
          cellOut = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(amountO);
          cellOut = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(amountN);
          cellOut = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(amountCl);
          cellOut = outRow.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
          cellOut.setCellValue(parser.calculateTheoreticalMass(formulaString, false));
        }
        rowCountOut++;
      }
      resultWorkbook.write(out);
      out.close();
      myxls.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void translateMzMineFilesToExcel(){
    double mzGrouping = 0.002d;
    double intCutoff = 1000000000d;
    double protonMass = 1.007276387d;
    String dirPath = "E:\\Marlene\\mzMinePeaks";
    
    File mzMineDir = new File(dirPath);
    File[] files = mzMineDir.listFiles();
    List<Double> mzValues = new ArrayList<Double>();
    for (int i=0; i!=files.length;i++){
      String path = files[i].getAbsolutePath();
      if (path.endsWith(".csv")){
        Vector<Double> foundMzs = this.parseMzMineExport(path);
        for (Double foundMz : foundMzs){
          Double existingMz = existsDoubleInList(foundMz,mzValues,mzGrouping);        
          if (existingMz>0){
            mzValues.remove(existingMz);
            double mean = (existingMz+foundMz)/2;
            mzValues.add(mean);
          }else
            mzValues.add(foundMz);
        }
      }
    }
    List<Double> candidateMzValues = new ArrayList<Double>();
    for (Double value : mzValues){
      double candidateSingle = value+protonMass;
      double candidateDouble = value*2+2*protonMass;
      Double existingSingle = existsDoubleInList(candidateSingle,candidateMzValues,mzGrouping);
      Double existingDouble = existsDoubleInList(candidateDouble,candidateMzValues,mzGrouping);
      if (existingSingle>0){
        candidateMzValues.remove(existingSingle);
        double mean = (existingSingle+candidateSingle)/2;
        candidateMzValues.add(mean);
      }else
        candidateMzValues.add(candidateSingle);
      if (existingDouble>0){
        candidateMzValues.remove(existingDouble);
        double mean = (existingDouble+candidateDouble)/2;
        candidateMzValues.add(mean);
      }else
        candidateMzValues.add(candidateDouble);
    }
    Collections.sort(candidateMzValues);
    // Juergen: now I parse known analytes and build the mean of C H O of the 7 nearest analytes - I will try to avoid duplicates
    // if the mass is higher than the highest mass - I select compounds with half of the mass - and duplicate the amounts
    try{
      InputStream myxls = new FileInputStream("E:\\Marlene\\Clematis_compounds.xls");
      HSSFWorkbook workbook     = new HSSFWorkbook(myxls);
      HSSFSheet sheet = workbook.getSheetAt(0);
      Hashtable<String,String> formulaHash = new Hashtable<String,String>();
      List<MassFormulaVO> massFormulas = new ArrayList<MassFormulaVO>();
      for (int rowCount=4;rowCount!=(256);rowCount++){
        HSSFRow row = sheet.getRow(rowCount); 
        String formulaString = "";
        String cString = "C";
        int amountC = 0;
        HSSFCell cell = row.getCell(1);
        Double value = getNumericCellValue(cell);
        if (value!=null) amountC = (int)Math.round(value);
        if (amountC > 0)
          formulaString += cString+String.valueOf(amountC);
        String hString = "H";
        int amountH = 0;
        cell = row.getCell(2);
        value = getNumericCellValue(cell);
        if (value!=null) amountH = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+hString+String.valueOf(amountH);
        String oString = "O";
        int amountO = 0;
        cell = row.getCell(3);
        value = getNumericCellValue(cell);
        if (value!=null) amountO = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+oString+String.valueOf(amountO);
        if (!formulaHash.containsKey(formulaString)){
          cell = row.getCell(5);
          double mass = getNumericCellValue(cell);
          formulaHash.put(formulaString, formulaString);
          MassFormulaVO formVO = new MassFormulaVO(mass,formulaString);
          massFormulas.add(formVO);
        }
      }
      Collections.sort(massFormulas, new GeneralComparator("at.tugraz.genome.TestClass$MassFormulaVO", "getMass", "java.lang.Double"));
      myxls.close();
      
      //now I have my formulasVOs - I will find now the 7 closest ones
      String quantFile = "E:\\Marlene\\201403_untargeted\\Clemantis_compounds_mzMine.xls";
      
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(quantFile));
      HSSFWorkbook resultWorkbook = new HSSFWorkbook();
      HSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
      HSSFSheet resultSheet = resultWorkbook.createSheet("Clemantis");
      int rowCountOut = 3;
      HSSFRow outRow = resultSheet.createRow(rowCountOut);
      HSSFCell label = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("Name");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("C");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("H");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("O");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("form[-H] name[-H]");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("tR (min)");
      label.setCellStyle(headerStyle);
      rowCountOut++;

      for (Double mz : candidateMzValues){
        outRow = resultSheet.createRow(rowCountOut);
        boolean isHalfValue = false;
        Double mzToSearch = new Double(mz);
        if (mz>2000d){
          mzToSearch = mz/2d;
          isHalfValue = true;
        }
        int position = 0;
        for (position = 0 ; position!=massFormulas.size(); position++){
          if (mzToSearch<massFormulas.get(position).getMass()) break;
        }
        // I have now the first position that is higher than the m/z value of interest - now have to get the 7 closest values
        Vector<MassFormulaVO> sevenClosestValues = new Vector<MassFormulaVO>();
        for (int i=position; (i!=(position+7)&&i!=massFormulas.size());i++){
          sevenClosestValues.add(massFormulas.get(i));
        }
        for (int i=(position-1); (i!=(position-8)&&i!=-1);i--){
          // first I have to find out if the Vector is already full
          // if it is full - add the element
          if (sevenClosestValues.size()<7) sevenClosestValues.add(massFormulas.get(i));
          else{
            // now, I have to find out if there is any value that has a higher difference to my mzValue
            // if yes, I have to throw out the value with the highest difference, if no - I can stop the for clause
            int farestPosition = -1;
            double maxDist = Math.abs(massFormulas.get(i).getMass()-mz);
            for (int j=0; j!=sevenClosestValues.size(); j++){
              double dist = Math.abs(sevenClosestValues.get(j).getMass()-mz);
              if (dist>maxDist){
                maxDist = dist;
                farestPosition = j;
              }
            }
            // there is one value in the Vector that is farer away
            if (farestPosition>-1){
              sevenClosestValues.remove(farestPosition);
              sevenClosestValues.add(massFormulas.get(i));
            }else{
              break;
            }
          }
        }
        //now calculate the mean of the elemental composition
        int amountC = 0;
        int amountH = 0;
        int amountO = 0;
        for (MassFormulaVO formVO : sevenClosestValues){
          StringTokenizer tokenizer = new StringTokenizer(formVO.getFormula()," ");
          amountC += Integer.parseInt(tokenizer.nextToken().substring(1));
          amountH += Integer.parseInt(tokenizer.nextToken().substring(1));
          amountO += Integer.parseInt(tokenizer.nextToken().substring(1));
        }
        if (isHalfValue){
          amountC = (int)Math.round((double)amountC*1.8d);
          amountH = amountH*2;
          amountO = (int)Math.round((double)amountO*2.22d);
        }
        amountC = amountC/sevenClosestValues.size();
        amountH = amountH/sevenClosestValues.size();
        amountO = amountO/sevenClosestValues.size();
        
        HSSFCell cellOut = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
        cellOut.setCellValue(Calculator.FormatNumber(mz, 3));
        cellOut = outRow.createCell(1,HSSFCell.CELL_TYPE_NUMERIC);
        cellOut.setCellValue(amountC);
        cellOut = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
        cellOut.setCellValue(amountH);
        cellOut = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
        cellOut.setCellValue(amountO);
        cellOut = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
        cellOut.setCellValue(mz);

        rowCountOut++;

      }
      resultWorkbook.write(out);
      out.close();
    }catch (Exception ex){
      ex.printStackTrace();
    }

  }
  
  
  public class MassFormulaVO {
    
    private String formula;
    private double mass;
    private String name;
    
    public MassFormulaVO(double mass, String formula){
      this.mass = mass;
      this.formula = formula;
    }

    public String getFormula()
    {
      return formula;
    }

    public double getMass()
    {
      return mass;
    }

    public void setName(String name)
    {
      this.name = name;
    }

    public String getName()
    {
      return name;
    }
    
    
  }
  
  private Double existsDoubleInList(Double foundMz, List<Double> list, double tol){
    Double existing = -1d;
    for (Double mz : list){
      if ((mz-tol)<foundMz && foundMz<(mz+tol)){
        existing = mz;
        break;
      }
    }
    return existing;
  }
  
  private Vector<Double> parseMzMineExport(String filePath){
    Vector<Double> masses = new Vector<Double>();
    try{
      LineNumberReader reader = new LineNumberReader(new FileReader(filePath));
      String line = reader.readLine();
      double maxArea = 0;
      while ((line = reader.readLine()) != null) {
        StringTokenizer tokenizer = new StringTokenizer(line,",");
        int maxTokens = tokenizer.countTokens();
        String mz = "";
        double mzNumeric = -1d;
        double area = -1;
        for (int i=0; i!=maxTokens;i++){
          String token = tokenizer.nextToken();
          if (i==1){
            mz = new String(token);
            mzNumeric = Double.parseDouble(mz);
          }
          if (i==(maxTokens-1)){
            area = Double.parseDouble(token);
          }
        }
        if (area>maxArea)maxArea = area;
//        if (area>cutoff){
        if (area>(maxArea/10d)){
          masses.add(mzNumeric);
        }
      }
    } catch (Exception ex){
      ex.printStackTrace();
    }
    return masses;
  }

  
  private void mergeNamedAndUntargetedExcel(){
    double mzGrouping = 0.002d;
    try{
      InputStream myxls = new FileInputStream("E:\\Marlene\\Clematis_compounds.xls");
      HSSFWorkbook workbook     = new HSSFWorkbook(myxls);
      HSSFSheet sheet = workbook.getSheetAt(0);
      List<MassFormulaVO> massFormulasNamed = new ArrayList<MassFormulaVO>();
      for (int rowCount=4;rowCount!=(256);rowCount++){
        HSSFRow row = sheet.getRow(rowCount); 
        HSSFCell cell = row.getCell(0);
        String name = cell.getStringCellValue();
        String formulaString = "";
        String cString = "C";
        int amountC = 0;
        cell = row.getCell(1);
        Double value = getNumericCellValue(cell);
        if (value!=null) amountC = (int)Math.round(value);
        if (amountC > 0)
          formulaString += cString+String.valueOf(amountC);
        String hString = "H";
        int amountH = 0;
        cell = row.getCell(2);
        value = getNumericCellValue(cell);
        if (value!=null) amountH = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+hString+String.valueOf(amountH);
        String oString = "O";
        int amountO = 0;
        cell = row.getCell(3);
        value = getNumericCellValue(cell);
        if (value!=null) amountO = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+oString+String.valueOf(amountO);
        String nString = "N";
        int amountN = 0;
        cell = row.getCell(4);
        value = getNumericCellValue(cell);
        if (value!=null) amountN = (int)Math.round(value);
        if (amountN > 0)
          formulaString += " "+nString+String.valueOf(amountN);
/*        String clString = "Cl";
        int amountCl = 0;
        cell = row.getCell(5);
        value = getNumericCellValue(cell);
        if (value!=null) amountCl = (int)Math.round(value);
        if (amountCl > 0)
          formulaString += " "+clString+String.valueOf(amountCl);*/
        cell = row.getCell(5);
        double mass = getNumericCellValue(cell);
        MassFormulaVO formVO = new MassFormulaVO(mass,formulaString);
        formVO.setName(name);
        massFormulasNamed.add(formVO);
      }
      Collections.sort(massFormulasNamed, new GeneralComparator("at.tugraz.genome.TestClass$MassFormulaVO", "getMass", "java.lang.Double"));
      myxls.close();
      
      List<MassFormulaVO> massFormulasUntargeted = new ArrayList<MassFormulaVO>();
      myxls = new FileInputStream("E:\\Marlene\\201403_untargeted\\Clemantis_compounds_mzMine.xls");
      workbook     = new HSSFWorkbook(myxls);
      sheet = workbook.getSheetAt(0);
      for (int rowCount=4;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
        HSSFRow row = sheet.getRow(rowCount); 
        HSSFCell cell = row.getCell(0);
        String name = String.valueOf(cell.getNumericCellValue());
        String formulaString = "";
        String cString = "C";
        int amountC = 0;
        cell = row.getCell(1);
        Double value = getNumericCellValue(cell);
        if (value!=null) amountC = (int)Math.round(value);
        if (amountC > 0)
          formulaString += cString+String.valueOf(amountC);
        String hString = "H";
        int amountH = 0;
        cell = row.getCell(2);
        value = getNumericCellValue(cell);
        if (value!=null) amountH = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+hString+String.valueOf(amountH);
        String oString = "O";
        int amountO = 0;
        cell = row.getCell(3);
        value = getNumericCellValue(cell);
        if (value!=null) amountO = (int)Math.round(value);
        if (amountH > 0)
          formulaString += " "+oString+String.valueOf(amountO);
        cell = row.getCell(4);
        double mass = getNumericCellValue(cell);
        MassFormulaVO formVO = new MassFormulaVO(mass,formulaString);
        formVO.setName(name);
        massFormulasUntargeted.add(formVO);
      }
      myxls.close();
      
      List<MassFormulaVO> finalList = new ArrayList<MassFormulaVO>();
      int currentNamed = 0;
      MassFormulaVO namedVO = massFormulasNamed.get(currentNamed);
      for (MassFormulaVO untVO : massFormulasUntargeted){
        boolean foundNamed = false;
        while (currentNamed<massFormulasNamed.size() && (namedVO.getMass()-mzGrouping)<untVO.getMass() && untVO.getMass()<(namedVO.getMass()+mzGrouping)){
          foundNamed = true;
          finalList.add(namedVO);
          currentNamed++;
          if (currentNamed>=massFormulasNamed.size()) break;
          namedVO = massFormulasNamed.get(currentNamed);
        }
        if (foundNamed){
          System.out.println("removing: "+untVO.getName());
          continue;
        }
        if (untVO.getMass()<namedVO.getMass()){
          finalList.add(untVO);
        }else{
          while (currentNamed<massFormulasNamed.size() && untVO.getMass()>=namedVO.getMass()){
            if ((namedVO.getMass()-mzGrouping)<untVO.getMass() && untVO.getMass()<(namedVO.getMass()+mzGrouping)) foundNamed = true;
            finalList.add(namedVO);
            currentNamed++;
            if (currentNamed>=massFormulasNamed.size()) break;
            namedVO = massFormulasNamed.get(currentNamed);            
          }
          while (currentNamed<massFormulasNamed.size() && (namedVO.getMass()-mzGrouping)<untVO.getMass() && untVO.getMass()<(namedVO.getMass()+mzGrouping)){
            foundNamed = true;
            finalList.add(namedVO);
            currentNamed++;
            if (currentNamed>=massFormulasNamed.size()) break;
            namedVO = massFormulasNamed.get(currentNamed);
          }
          if (!foundNamed) finalList.add(untVO);
          else System.out.println("removing: "+untVO.getName());
        }
      }
      
      String quantFile = "E:\\Marlene\\201403_untargeted\\Clemantis_compounds.xls";
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(quantFile));
      HSSFWorkbook resultWorkbook = new HSSFWorkbook();
      HSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
      HSSFSheet resultSheet = resultWorkbook.createSheet("Clemantis");
      int rowCountOut = 3;
      HSSFRow outRow = resultSheet.createRow(rowCountOut);
      HSSFCell label = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("Name");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("C");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("H");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("O");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("N");
      label.setCellStyle(headerStyle);
//      label = outRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
//      label.setCellValue("Cl");
//      label.setCellStyle(headerStyle);

      label = outRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("form[-H] name[-H]");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("tR (min)");
      label.setCellStyle(headerStyle);
      rowCountOut++;
      System.out.println("FinalListSize: "+finalList.size());
      for (MassFormulaVO formVO : finalList){
        StringTokenizer tokenizer = new StringTokenizer(formVO.getFormula()," ");
        int amountC = 0;
        int amountH = 0;
        int amountO = 0;
        int amountN = 0;
//        int amountCl = 0;
        while (tokenizer.hasMoreTokens()){
          String token  = tokenizer.nextToken();
//          if (token.startsWith("Cl")) amountCl = Integer.parseInt(token.substring(2));
/*          else */if (token.startsWith("N")) amountN = Integer.parseInt(token.substring(1));
          else if (token.startsWith("O")) amountO = Integer.parseInt(token.substring(1));
          else if (token.startsWith("H")) amountH = Integer.parseInt(token.substring(1));
          else if (token.startsWith("C")) amountC = Integer.parseInt(token.substring(1));
        }
        outRow = resultSheet.createRow(rowCountOut);
        HSSFCell cellOut = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
        cellOut.setCellValue(formVO.getName());
        cellOut = outRow.createCell(1,HSSFCell.CELL_TYPE_NUMERIC);
        cellOut.setCellValue(amountC);
        cellOut = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
        cellOut.setCellValue(amountH);
        cellOut = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
        cellOut.setCellValue(amountO);
        cellOut = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
        cellOut.setCellValue(amountN);
//        cellOut = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
//        cellOut.setCellValue(amountCl);
        
        cellOut = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
        cellOut.setCellValue(formVO.getMass());

        rowCountOut++;
      }      
      resultWorkbook.write(out);
      out.close();
    }catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void detectMsn(){
    try {
      /*      LipidParameterSet param = new LipidParameterSet(885.549793749f, "38", 4, "H", "", "C47 H83 O13 P1", "-H",1);

      CgProbe probe2 = new CgProbe(0,1);
      probe2.AreaStatus = CgAreaStatus.OK;
      probe2.Area = 16958208f;
      probe2.AreaError = 2232620.25f;
      probe2.Background = 5865.395996f;
      probe2.Peak = 567.4219971f;
      probe2.LowerValley = 534.2869873f;
      probe2.UpperValley = 578.2659912f;
      probe2.Mz = 885.5488281f;
      probe2.LowerMzBand = 0.013f;
      probe2.UpperMzBand = 0.013f;
      probe2.isotopeNumber = 0;
      param.AddProbe(probe2);
      
      CgProbe probe1 = new CgProbe(0,1);
      probe1.AreaStatus = CgAreaStatus.OK;
      probe1.Area = 259979168f;
      probe1.AreaError = 36839124f;
      probe1.Background = 174605.2969f;
      probe1.Peak = 619.8410034f;
      probe1.LowerValley = 600.0460205f;
      probe1.UpperValley = 647.0319824f;
      probe1.Mz = 885.5488281f;
      probe1.LowerMzBand = 0.013f;
      probe1.UpperMzBand = 0.013f;
      probe1.isotopeNumber = 0;
      param.AddProbe(probe1);*/

//      LipidParameterSet param = new LipidParameterSet(852.80143449f, "50", 0, "NH4", "", "C53 H102 O6", "+NH4",1);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 2458533376f;
//      probe1.AreaError = 171287568f;
//      probe1.Background = 1073633.875f;
//      probe1.Peak = 2509.98999f;
//      probe1.LowerValley = 2486.98999f;
//      probe1.UpperValley = 2541.439941f;
//      probe1.Mz = 852.805542f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);

//      LipidParameterSet param = new LipidParameterSet(768.707534286f, "44", 0, "NH4", "", "C47 H90 O6", "+NH4",1);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 101906232f;
//      probe1.AreaError = 5509909.5f;
//      probe1.Background = 41727.40234f;
//      probe1.Peak = 2321.86010742187f;
//      probe1.LowerValley = 2302.55004882812f;
//      probe1.UpperValley = 2356.080078125f;
//      probe1.Mz = 768.709594726562f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);

//      LipidParameterSet param = new LipidParameterSet(850.7857666f, "50", 1, "NH4", "", "C57 H106 O6", "+NH4",1);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 2728108.125f;
//      probe1.AreaError = 145217.828125f;
//      probe1.Background = 8609.3466796875f;
//      probe1.Peak = 2458.25f;
//      probe1.LowerValley = 2434.59008789062f;
//      probe1.UpperValley = 2489.4599609375f;
//      probe1.Mz = 850.787841796875f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
      
//      LipidParameterSet param = new LipidParameterSet(878.81708f, "52", 1, "NH4", "", "C45 H86 O6", "+NH4",1);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 93738096f;
//      probe1.AreaError = 5153103f;
//      probe1.Background = 42592.51172f;
//      probe1.Peak = 2515.310059f;
//      probe1.LowerValley = 2491.949951f;
//      probe1.UpperValley = 2548.459961f;
//      probe1.Mz = 878.8201294f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
      
      
//      LipidParameterSet param = new LipidParameterSet(885.5498046875f, "38", 4, "H", "", "H83 P1 O13 N1 C47", "-H",1);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 93738096f;
//      probe1.AreaError = 5153103f;
//      probe1.Background = 42592.51172f;
//      probe1.Peak = 1555.13000488281f;
//      probe1.LowerValley = 1538.72998046875f;
//      probe1.UpperValley = 1576.7099609375f;
//      probe1.Mz = 885.548828125f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);

      
//      LipidParameterSet param = new LipidParameterSet(608.524841308593f, "36", 3, "-H", "", "H68 O5 C35", "-H",1);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 93738096f;
//      probe1.AreaError = 5153103f;
//      probe1.Background = 42592.51172f;
//      probe1.Peak = 1719.85998535156f;
//      probe1.LowerValley = 1690.17004394531f;
//      probe1.UpperValley = 1787.47998046875f;
//      probe1.Mz = 608.522888183593f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
      
      LipidParameterSet param = new LipidParameterSet(874.785827636718f, "52", 3, 0, "HH4", "", "H100 O6 C55", "+NH4",1);

      CgProbe probe1 = new CgProbe(0,1);
      probe1.AreaStatus = CgAreaStatus.OK;
      probe1.Area = 3016949760f;
      probe1.AreaError = 743906944f;
      probe1.Background = 1286404.5f;
      probe1.Peak = 1277.55004882812f;
      probe1.LowerValley = 1262.55004882812f;
      probe1.UpperValley = 1292.75f;
      probe1.Mz = 874.789855957031f;
      probe1.LowerMzBand = 0.013f;
      probe1.UpperMzBand = 0.013f;
      probe1.isotopeNumber = 0;
      param.AddProbe(probe1);

      
      //String[] chromPaths = StringUtils.getChromFilePaths("E:\\lipidomicsMS2\\20140220\\001_LD1.chrom");
      //String[] chromPaths = StringUtils.getChromFilePaths("E:\\lipidomicsMS2\\20140220\\002_LD2.chrom");
      //String[] chromPaths = StringUtils.getChromFilePaths("D:\\Experiment1\\LipidBLAST\\002_liver_pos_CID50.chrom");
      ////String[] chromPaths = StringUtils.getChromFilePaths("E:\\lipidomicsMS2\\20141030\\neg\\004_neg_CID_50.chrom");
      String[] chromPaths = StringUtils.getChromFilePaths("D:\\Evelyn\\20171006\\SRM1950_1.chrom");
      
      System.out.println(chromPaths[1]);
      LipidomicsAnalyzer lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
      setStandardParameters(lAnalyzer);
      MSnAnalyzer analyzer = new MSnAnalyzer(null,"TG","NH4",param,lAnalyzer,null,false,true,true);
//      MSnDebugVO debugInfo = analyzer.getDebugInfo();
//      System.out.println("Debug-Status: "+debugInfo.getStatus());
//      Hashtable<String, Integer> discHeads = debugInfo.getDiscardedHeadGroupFragments();
//      System.out.println("Discarded-Heads-Size: "+discHeads.size());
//      for (String fragment: discHeads.keySet()){
//        System.out.println("Head-Discarded: "+fragment+"\t"+discHeads.get(fragment));
//      }
//      Hashtable<String,IntensityRuleVO> violHeadRules = debugInfo.getViolatedHeadRules();
//      System.out.println("Violated head rules size: "+violHeadRules.size());
//      for (IntensityRuleVO ruleVO : violHeadRules.values()){
//        System.out.println(ruleVO.getReadableRuleInterpretation());
//      }
//      Hashtable<String,Hashtable<String,Object>> violChainRules = debugInfo.getViolatedChainRules();
//      System.out.println("Violated chain rules size: "+violChainRules.size());
//      for (String faName : violChainRules.keySet()){
//        Hashtable<String,Object> violRule = violChainRules.get(faName);
//        for (String ruleName : violRule.keySet()){
//          Object rule = violRule.get(ruleName);
//          if (rule instanceof Integer){
//            int status = (Integer)rule;
//            System.out.println("Chain Discarded: "+faName+"\t"+ruleName+"\t"+status);
//          } else if (rule instanceof IntensityRuleVO){
//            IntensityRuleVO ruleVO = (IntensityRuleVO)rule;
//            System.out.println("Violated chain rule: "+faName+"\t"+ruleVO.getReadableRuleInterpretation());
//          }
//        }
//      }
//      Hashtable<String,Integer> violCombis = debugInfo.getViolatedCombinations();
//      System.out.println("Violated combinations: "+violCombis.size());
//      for (String combi : violCombis.keySet()){
//        System.out.println("Violated combi: "+combi+"\t"+violCombis.get(combi));
//      }
//      System.out.println("Spectrum sufficiently covered: "+debugInfo.isSpectrumCoverageFulfilled());
//      Hashtable<String,Hashtable<String,IntensityRuleVO>> unfulfilledPosRules = debugInfo.getUnfulfilledPositionRules();
//      System.out.println("Unfulfilled rules: "+unfulfilledPosRules.size());
//      for (String combiName : unfulfilledPosRules.keySet()){
//        Hashtable<String,IntensityRuleVO> unfulfilled = unfulfilledPosRules.get(combiName);
//        for (IntensityRuleVO ruleVO : unfulfilled.values()){
//          System.out.println(combiName+":\tNOT "+ruleVO.getReadableRuleInterpretation());
//        }
//      }
//      Hashtable<String,Vector<Vector<IntensityRuleVO>>> contradictingPositionRules  = debugInfo.getContradictingPositionRules();
//      System.out.println("Any contradicting position rules: "+contradictingPositionRules.size());
//      for (String combiName : contradictingPositionRules.keySet()){
//        Vector<Vector<IntensityRuleVO>> rulePairs = contradictingPositionRules.get(combiName);
//        for (Vector<IntensityRuleVO> rulePair : rulePairs){
//          IntensityRuleVO rule1 = rulePair.get(0);
//          IntensityRuleVO rule2 = rulePair.get(1);
//          System.out.println(combiName+":\t"+rule1.getReadableRuleInterpretation()+"\t!=\t"+rule2.getReadableRuleInterpretation());
//        }
//      }
      
      
      System.out.println("Status: "+analyzer.checkStatus());
      if (analyzer.checkStatus()!=LipidomicsMSnSet.NO_MSN_PRESENT && analyzer.checkStatus()!=LipidomicsMSnSet.DISCARD_HIT){
        Hashtable<String,CgProbe> headFragments = analyzer.getHeadGroupFragments();
        for (String fragmentName : headFragments.keySet()){
          CgProbe probe = headFragments.get(fragmentName);
          System.out.println(fragmentName+": \t"+probe.Mz+"\t"+probe.Area);
        }
        Hashtable<String,Hashtable<String,CgProbe>> chainFragments = analyzer.getChainFragments();
        for (String chain : chainFragments.keySet()){
          System.out.println("1. "+chain);
          Hashtable<String,CgProbe> frags = chainFragments.get(chain);
          for (String frag : frags.keySet()){
            CgProbe probe = frags.get(frag);
            System.out.println(frag+": \t"+probe.Mz+"\t"+probe.Area);
          }
        }
//        Hashtable<String,Double> intensities = analyzer.getChainCombinationRelativeAreas();
//        for (String key : intensities.keySet()){
//          System.out.println(key+" "+intensities.get(key));
//        }
      }
      printResults(analyzer.getResult());
      
//      param = new LipidParameterSet(878.81708f, "52", 1, "NH4", "", "C45 H86 O6", "+NH4",1);
//      probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 93738096f;
//      probe1.AreaError = 5153103f;
//      probe1.Background = 42592.51172f;
//      probe1.Peak = 2515.310059f;
//      probe1.LowerValley = 2491.949951f;
//      probe1.UpperValley = 2548.459961f;
//      probe1.Mz = 878.8201294f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
//      chromPaths = StringUtils.getChromFilePaths("E:\\lipidomicsMS2\\20140220\\001_LD1.chrom");
//      System.out.println(chromPaths[1]);
//      lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0]);
//      setStandardParameters(lAnalyzer);
//      analyzer = new MSnAnalyzer(null,"TG","NH4",param,lAnalyzer,false);
//      System.out.println("Status: "+analyzer.checkStatus());
//      printResults(analyzer.getResult());
//      System.out.println("------------------------------------");
//      
////   -------------------------------------------------------
//      
//      param = new LipidParameterSet(861.549793749f, "36", 2, "H", "", "C45 H83 O13 P1", "-H",1);
//
//      probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 93738096f;
//      probe1.AreaError = 5153103f;
//      probe1.Background = 42592.51172f;
//      probe1.Peak = 1422.73999023437f;
//      probe1.LowerValley = 1357.59997558593f;
//      probe1.UpperValley = 1455.08996582031f;
//      probe1.Mz = 861.5478515625f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);      
//      chromPaths = StringUtils.getChromFilePaths("E:\\lipidomicsMS2\\20140218\\neg\\030_Liver1.chrom");
//      lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0]);
//      setStandardParameters(lAnalyzer);
//      analyzer = new MSnAnalyzer("tmp","PI","H",param,lAnalyzer,false);
//      System.out.println("Status: "+analyzer.checkStatus());
//      printResults(analyzer.getResult());
//      System.out.println("------------------------------------");
//
//      param = new LipidParameterSet(878.81708f, "52", 1, "NH4", "", "C45 H86 O6", "+NH4",1);
//      probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 93738096f;
//      probe1.AreaError = 5153103f;
//      probe1.Background = 42592.51172f;
//      probe1.Peak = 2515.310059f;
//      probe1.LowerValley = 2491.949951f;
//      probe1.UpperValley = 2548.459961f;
//      probe1.Mz = 878.8201294f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
//      chromPaths = StringUtils.getChromFilePaths("E:\\lipidomicsMS2\\20140220\\001_LD1.chrom");
//      System.out.println(chromPaths[1]);
//      lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0]);
//      setStandardParameters(lAnalyzer);
//      try{
//        analyzer = new MSnAnalyzer("tmp","TG","NH4",param,lAnalyzer,false);
//        System.out.println("Status: "+analyzer.checkStatus());
//        printResults(analyzer.getResult());
//      } catch (Exception ex){
//        System.out.println("In the tmp folder, there are no rules for TG: "+ex.getMessage());
//      }
//      System.out.println("------------------------------------");
//
//      param = new LipidParameterSet(861.549793749f, "36", 2, "H", "", "C45 H83 O13 P1", "-H",1);
//
//      probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 93738096f;
//      probe1.AreaError = 5153103f;
//      probe1.Background = 42592.51172f;
//      probe1.Peak = 1422.73999023437f;
//      probe1.LowerValley = 1357.59997558593f;
//      probe1.UpperValley = 1455.08996582031f;
//      probe1.Mz = 861.5478515625f;
//      probe1.LowerMzBand = 0.013f;
//      probe1.UpperMzBand = 0.013f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);      
//      chromPaths = StringUtils.getChromFilePaths("E:\\lipidomicsMS2\\20140218\\neg\\030_Liver1.chrom");
//      lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0]);
//      setStandardParameters(lAnalyzer);
//      analyzer = new MSnAnalyzer(null,"PI","H",param,lAnalyzer,false);
//      System.out.println("Status: "+analyzer.checkStatus());
//      printResults(analyzer.getResult());
//
//
////      Hashtable<Integer,String> definedPositions = new Hashtable<Integer,String>();
////      definedPositions.put(1, "18:0");
////      Vector<String> unassignedFAs = new Vector<String>();
////      unassignedFAs.add("17:0");
////      unassignedFAs.add("19:0");
////      unassignedFAs.add("20:0");
////      Vector<String> names = MSnAnalyzer.getPermutedChainPositionNames(definedPositions, unassignedFAs);
////      System.out.println(names);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  private void detectPlasmalogens(){
/*    LipidParameterSet param = new LipidParameterSet(728.559416177f, "36", 2, "H", "", "C41 H80 O7 P1 N1", "-H",1);
    CgProbe probe1 = new CgProbe(0,1);
    probe1.AreaStatus = CgAreaStatus.OK;
    probe1.Area = 383439840f;
    probe1.AreaError = 4248196.5f;
    probe1.Background = 404001.09375f;
    probe1.Peak = 843.255981445312f;
    probe1.LowerValley = 831.366027832031f;
    probe1.UpperValley = 852.598022460937f;
    probe1.Mz = 728.559387207031f;
    probe1.LowerMzBand = 0.015f;
    probe1.UpperMzBand = 0.015f;
    probe1.isotopeNumber = 0;*/
    
    LipidParameterSet param = new LipidParameterSet(728.559387207031f, "36", 2, 0, "-H", "", "C41 H80 O7 P1 N1", "-H",1);
    CgProbe probe1 = new CgProbe(0,1);
    probe1.AreaStatus = CgAreaStatus.OK;
    probe1.Area = 383439840f;
    probe1.AreaError = 4248196.5f;
    probe1.Background = 404001.09375f;
    probe1.Peak = 843.255981445312f;
    probe1.LowerValley = 831.366027832031f;
    probe1.UpperValley = 852.598022460937f;
    probe1.Mz = 728.559387207031f;
    probe1.LowerMzBand = 0.015f;
    probe1.UpperMzBand = 0.015f;
    probe1.isotopeNumber = 0;

    try {
      param.AddProbe(probe1);
      String[] chromPaths = StringUtils.getChromFilePaths("E:\\Robin\\20140819\\Control\\10uL_21Mar14_L4440_Negative_MS2_32-97-32_300uLmin_01.chrom");
      LipidomicsAnalyzer lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
      setStandardParameters(lAnalyzer);
      MSnAnalyzer analyzer = new MSnAnalyzer(null,"O-PE","-H",param,lAnalyzer,null,false,true,false);
      System.out.println("Status: "+analyzer.checkStatus());
      printResults(analyzer.getResult());
      LipidomicsMSnSet set = (LipidomicsMSnSet)analyzer.getResult();
      for (String chainCombi : set.getChainIntensityRules().keySet()) System.out.println(chainCombi);
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }      

  }
  
  private void alkylAlkenylCombinations(){
//    Vector<Vector<Integer>> combis = FragmentCalculator.alkylAlkenylPermutations(3, 2, 1);
//    for (Vector<Integer> combi : combis) System.out.println(combi);
  }
  
  private void printResults(LipidParameterSet result) throws LipidCombinameEncodingException{
    if (result instanceof LipidomicsMSnSet){
      LipidomicsMSnSet resultMSn = (LipidomicsMSnSet)result;
      System.out.println("Head-group size: "+resultMSn.getHeadGroupFragments().size());
      for (Object nameObject : resultMSn.getMSnIdentificationNames()){
        if (nameObject instanceof Vector){
          for (String name : (Vector<String>)nameObject){
            System.out.println("Proposed composition (Vec): "+name);              
          }
          System.out.println("Rel Area: "+resultMSn.getRelativeIntensity(((Vector<String>)nameObject).get(0)));
        }else{
          System.out.println("Proposed composition: "+(String)nameObject+" ; "+resultMSn.getRelativeIntensity((String)nameObject));
        }
      }
    }else{
      System.out.println(result.getNameStringWithoutRt());
    }    
  }
  
  private void setStandardParameters(LipidomicsAnalyzer lAnalyzer){
    lAnalyzer.set3DParameters(LipidomicsConstants.getChromSmoothRange(),LipidomicsConstants.getChromSmoothRepeats(),
        LipidomicsConstants.removeIfOtherIsotopePresent(),LipidomicsConstants.useNoiseCutoff(), LipidomicsConstants.getNoiseCutoffDeviationValue(), 
        LipidomicsConstants.getMinimumRelativeIntensity(), LipidomicsConstants.getScanStep(),
        LipidomicsConstants.getProfileMzRange(), LipidomicsConstants.getProfileTimeTolerance_(),
        LipidomicsConstants.getProfileIntThreshold_(), LipidomicsConstants.getBroaderProfileTimeTolerance_(),
        LipidomicsConstants.getProfileSmoothRange(),LipidomicsConstants.getProfileSmoothRepeats(),
        LipidomicsConstants.getProfileMeanSmoothRepeats(), LipidomicsConstants.getProfileMzMinRange(),
        LipidomicsConstants.getProfileSteepnessChange1(),LipidomicsConstants.getProfileSteepnessChange2(), 
        LipidomicsConstants.getProfileIntensityCutoff1(), LipidomicsConstants.getProfileIntensityCutoff2(), 
        LipidomicsConstants.getProfileGeneralIntCutoff(),LipidomicsConstants.getProfilePeakAcceptanceRange(), 
        LipidomicsConstants.getProfileSmoothingCorrection(), LipidomicsConstants.getProfileMaxRange(),
        LipidomicsConstants.getSmallChromMzRange(),LipidomicsConstants.getSmallChromSmoothRepeats(),
        LipidomicsConstants.getSmallChromMeanSmoothRepeats(), LipidomicsConstants.getSmallChromSmoothRange(),
        LipidomicsConstants.getSmallChromIntensityCutoff(),LipidomicsConstants.getBroadChromSmoothRepeats(),
        LipidomicsConstants.getBroadChromMeanSmoothRepeats(), LipidomicsConstants.getBroadChromSmoothRange(),
        LipidomicsConstants.getBroadChromIntensityCutoff(), LipidomicsConstants.getBroadChromSteepnessChangeNoSmall(),
        LipidomicsConstants.getBroadIntensityCutoffNoSmall(),
        LipidomicsConstants.getFinalProbeTimeCompTolerance(), LipidomicsConstants.getFinalProbeMzCompTolerance(),
        LipidomicsConstants.getOverlapDistanceDeviationFactor(), LipidomicsConstants.getOverlapPossibleIntensityThreshold(), 
        LipidomicsConstants.getOverlapSureIntensityThreshold(), LipidomicsConstants.getOverlapPeakDistanceDivisor(),
        LipidomicsConstants.getOverlapFullDistanceDivisor(), LipidomicsConstants.getPeakDiscardingAreaFactor(),
        LipidomicsConstants.getIsotopeInBetweenTime(),LipidomicsConstants.getIsoInBetweenAreaFactor(),LipidomicsConstants.getIsoInBetweenMaxTimeDistance(),
        LipidomicsConstants.getIsoNearNormalProbeTime(),LipidomicsConstants.getRelativeAreaCutoff(),
        LipidomicsConstants.getRelativeFarAreaCutoff(),LipidomicsConstants.getRelativeFarAreaTimeSpace(),
        LipidomicsConstants.getRelativeIsoInBetweenCutoff(),LipidomicsConstants.getTwinPeakMzTolerance(),
        LipidomicsConstants.getClosePeakTimeTolerance(),LipidomicsConstants.getTwinInBetweenCutoff(), LipidomicsConstants.getUnionInBetweenCutoff(),
        LipidomicsConstants.getMs2MzTolerance());
  }
  
  private Metadata getMetadata() throws Exception {
    //MZTabDescription tabDescription = new MZTabDescription(MZTabDescription.Mode.Summary, MZTabDescription.Type.Identification);
    //tabDescription.setId("PRIDE_1234");
    //Metadata mtd = new Metadata(tabDescription);
    Metadata mtd = new Metadata();

    mtd.setTitle("My first test experiment");
    mtd.setDescription("An experiment investigating the effects of Il-6.");
    //the following lines were not updated in the course of the mzTab 1.1 update
    /*
    mtd.addSampleProcessingParam(1, new CVParam("SEP", "SEP:00173", "SDS PAGE", null));
    mtd.addSampleProcessingParam(2, new CVParam("SEP", "SEP:00142", "enzyme digestion", null));
    mtd.addSampleProcessingParam(2, new CVParam("MS", "MS:1001251", "Trypsin", null));

    mtd.addInstrumentName(1, new CVParam("MS", "MS:100049", "LTQ Orbitrap", null));
    mtd.addInstrumentName(2, new CVParam("MS", "MS:1000031", "Instrument model", "name of the instrument not included in the CV"));
    mtd.addInstrumentSource(1, new CVParam("MS", "MS:1000073", "ESI", null));
    mtd.addInstrumentSource(2, new CVParam("MS", "MS:1000598", "ETD", null));
    mtd.addInstrumentAnalyzer(1, new CVParam("MS", "MS:1000291", "linear ion trap", null));
    mtd.addInstrumentAnalyzer(2, new CVParam("MS", "MS:1000484", "orbitrap", null));
    mtd.addInstrumentDetector(1, new CVParam("MS", "MS:1000253", "electron multiplier", null));
    mtd.addInstrumentDetector(2, new CVParam("MS", "MS:1000348", "focal plane collector", null));

    mtd.addSoftwareParam(1, new CVParam("MS", "MS:1001207", "Mascot", "2.3"));
    mtd.addSoftwareParam(2, new CVParam("MS", "MS:1001561", "Scaffold", "1.0"));
    mtd.addSoftwareSetting(1, "Fragment tolerance = 0.1Da");
    mtd.addSoftwareSetting(1, "Parent tolerance = 0.5Da");

    mtd.addFalseDiscoveryRateParam(new CVParam("MS", "MS:1001364", "pep:global FDR", "0.01"));
    mtd.addFalseDiscoveryRateParam(new CVParam("MS", "MS:1001214", "pep:global FDR", "0.08"));

    mtd.addPublicationItem(1, PublicationItem.Type.PUBMED, "21063943");
    mtd.addPublicationItem(1, PublicationItem.Type.DOI, "10.1007/978-1-60761-987-1_6");
    mtd.addPublicationItem(2, PublicationItem.Type.PUBMED, "20615486");
    mtd.addPublicationItem(2, PublicationItem.Type.DOI, "10.1016/j.jprot.2010.06.008");

    mtd.addContactName(1, "James D. Watson");
    mtd.addContactName(2, "Francis Crick");
    mtd.addContactAffiliation(1, "Cambridge University, UK");
    mtd.addContactAffiliation(2, "Cambridge University, UK");
    mtd.addContactEmail(1, "watson@cam.ac.uk");
    mtd.addContactEmail(2, "crick@cam.ac.uk");

////    mtd.addUri(new URI("http://www.ebi.ac.uk/pride/url/to/experiment"));
////    mtd.addUri(new URI("http://proteomecentral.proteomexchange.org/cgi/GetDataset"));

    mtd.addFixedModParam(1, new CVParam("UNIMOD", "UNIMOD:4", "Carbamidomethyl", null));
    mtd.addFixedModSite(1, "M");
    mtd.addFixedModParam(2, new CVParam("UNIMOD", "UNIMOD:35", "Oxidation", null));
    mtd.addFixedModSite(2, "N-term");
    mtd.addFixedModParam(3, new CVParam("UNIMOD", "UNIMOD:1", "Acetyl", null));
    mtd.addFixedModPosition(3, "Protein C-term");

    mtd.addVariableModParam(1, new CVParam("UNIMOD", "UNIMOD:21", "Phospho", null));
    mtd.addVariableModSite(1, "M");
    mtd.addVariableModParam(2, new CVParam("UNIMOD", "UNIMOD:35", "Oxidation", null));
    mtd.addVariableModSite(2, "N-term");
    mtd.addVariableModParam(3, new CVParam("UNIMOD", "UNIMOD:1", "Acetyl", null));
    mtd.addVariableModPosition(3, "Protein C-term");

    mtd.setQuantificationMethod(new CVParam("MS", "MS:1001837", "iTRAQ quantitation analysis", null));
    mtd.setProteinQuantificationUnit(new CVParam("PRIDE", "PRIDE:0000395", "Ratio", null));
    mtd.setPeptideQuantificationUnit(new CVParam("PRIDE", "PRIDE:0000395", "Ratio", null));
    mtd.setSmallMoleculeQuantificationUnit(new CVParam("PRIDE", "PRIDE:0000395", "Ratio", null));

    mtd.addMsRunFormat(1, new CVParam("MS", "MS:1000584", "mzML file", null));
    mtd.addMsRunFormat(2, new CVParam("MS", "MS:1001062", "Mascot MGF file", null));
    mtd.addMsRunLocation(1, new URL("file://C:\\path\\to\\my\\file"));
    mtd.addMsRunLocation(2, new URL("ftp://ftp.ebi.ac.uk/path/to/file"));
    mtd.addMsRunIdFormat(1, new CVParam("MS", "MS:1001530", "mzML unique identifier", null));
    mtd.addMsRunFragmentationMethod(1, new CVParam("MS", "MS:1000133", "CID", null));
    mtd.addMsRunFragmentationMethod(2, new CVParam("MS", "MS:1000422", "HCD", null));

    mtd.addCustom(new UserParam("MS operator", "Florian"));

    mtd.addSampleSpecies(1, new CVParam("NEWT", "9606", "Homo sapiens (Human)", null));
    mtd.addSampleSpecies(1, new CVParam("NEWT", "573824", "Human rhinovirus 1", null));
    mtd.addSampleSpecies(2, new CVParam("NEWT", "9606", "Homo sapiens (Human)", null));
    mtd.addSampleSpecies(2, new CVParam("NEWT", "12130", "Human rhinovirus 2", null));
    mtd.addSampleTissue(1, new CVParam("BTO", "BTO:0000759", "liver", null));
    mtd.addSampleCellType(1, new CVParam("CL", "CL:0000182", "hepatocyte", null));
    mtd.addSampleDisease(1, new CVParam("DOID", "DOID:684", "hepatocellular carcinoma", null));
    mtd.addSampleDisease(1, new CVParam("DOID", "DOID:9451", "alcoholic fatty liver", null));
    mtd.addSampleDescription(1, "Hepatocellular carcinoma samples.");
    mtd.addSampleDescription(2, "Healthy control samples.");
    mtd.addSampleCustom(1, new UserParam("Extraction date", "2011-12-21"));
    mtd.addSampleCustom(1, new UserParam("Extraction reason", "liver biopsy"));

    Sample sample1 = mtd.getSampleMap().get(1);
    Sample sample2 = mtd.getSampleMap().get(2);
    mtd.addAssayQuantificationReagent(1, new CVParam("PRIDE", "PRIDE:0000114", "iTRAQ reagent", "114"));
    mtd.addAssayQuantificationReagent(2, new CVParam("PRIDE", "PRIDE:0000115", "iTRAQ reagent", "115"));
    mtd.addAssayQuantificationReagent(1, new CVParam("PRIDE", "MS:1002038", "unlabeled sample", null));
    mtd.addAssaySample(1, sample1);
    mtd.addAssaySample(2, sample2);

    mtd.addAssayQuantificationModParam(2, 1, new CVParam("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
    mtd.addAssayQuantificationModParam(2, 2, new CVParam("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
    mtd.addAssayQuantificationModSite(2, 1, "R");
    mtd.addAssayQuantificationModSite(2, 2, "K");
    mtd.addAssayQuantificationModPosition(2, 1, "Anywhere");
    mtd.addAssayQuantificationModPosition(2, 2, "Anywhere");

    MsRun msRun1 = mtd.getMsRunMap().get(1);
    mtd.addAssayMsRun(1, msRun1);

    Assay assay1 = mtd.getAssayMap().get(1);
    Assay assay2 = mtd.getAssayMap().get(2);
    mtd.addStudyVariableAssay(1, assay1);
    mtd.addStudyVariableAssay(1, assay2);

    mtd.addStudyVariableSample(1, sample1);
    mtd.addStudyVariableDescription(1, "description Group B (spike-in 0.74 fmol/uL)");

    mtd.addCVLabel(1, "MS");
    mtd.addCVFullName(1, "MS");
    mtd.addCVVersion(1, "3.54.0");
    mtd.addCVURL(1, "http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");

    mtd.addProteinColUnit(ProteinColumn.RELIABILITY, new CVParam("MS", "MS:00001231", "PeptideProphet:Score", null));

    MZTabColumnFactory peptideFactory = MZTabColumnFactory.getInstance(Section.Peptide);
    PeptideColumn peptideColumn = (PeptideColumn) peptideFactory.findColumnByHeader("retention_time");
    mtd.addPeptideColUnit(peptideColumn, new CVParam("UO", "UO:0000031", "minute", null));

    mtd.addPSMColUnit(PSMColumn.RETENTION_TIME, new CVParam("UO", "UO:0000031", "minute", null));
    mtd.addSmallMoleculeColUnit(SmallMoleculeColumn.RETENTION_TIME, new CVParam("UO", "UO:0000031", "minute", null));
*/
    return mtd;
}

//the following lines were not updated in the course of the mzTab 1.1 update
/*private MZTabColumnFactory getSMH(Metadata metadata) {
    MZTabColumnFactory factory = MZTabColumnFactory.getInstance(Section.Small_Molecule);

    // add optional columns which have stable order.
    factory.addReliabilityOptionalColumn();
    factory.addURIOptionalColumn();

    // add optional columns which have stable order.
    factory.addOptionalColumn(SmallMoleculeColumn.SEARCH_ENGINE_SCORE, metadata.getMsRunMap().get(1));

    // add abundance columns which locate the end of table.
    factory.addAbundanceOptionalColumn(metadata.getAssayMap().get(1));
    factory.addAbundanceOptionalColumn(metadata.getStudyVariableMap().get(1));
    factory.addAbundanceOptionalColumn(metadata.getAssayMap().get(2));

    // add user defined optional columns
    factory.addOptionalColumn(metadata.getMsRunMap().get(1), "my_value", String.class);
    CVParam param = new CVParam("MS", "MS:1002217", "decoy peptide", null);
    factory.addOptionalColumn(param, String.class);

    return factory;
}

private SmallMolecule getRecord(Metadata metadata, MZTabColumnFactory factory) {
    SmallMolecule sm = new SmallMolecule(factory, metadata);
    sm.setIdentifier("CID:00027395");
    sm.setChemicalFormula("C17H20N4O2");
    sm.setSmiles("C1=CC=C(C=C1)CCNC(=O)CCNNC(=O)C2=CC=NC=C2");
    sm.setInchiKey("QXBMEGUKVLFJAM-UHFFFAOYSA-N");
    sm.setDescription("N-(2-phenylethyl)-3-[2-(pyridine-4-carbonyl)hydrazinyl]propanamide");
    sm.setExpMassToCharge("1234.4");
    sm.setCalcMassToCharge("1234.5");
    sm.setCharge("2");
    sm.setRetentionTime("10.2|11.5");
    sm.setTaxid("10116");
    sm.setSpecies("Rattus norvegicus (Rat)");
    sm.setDatabase("UniProtKB");
    sm.setDatabaseVersion("2011_11");
    sm.setReliability("2");
    sm.setURI("http://www.ebi.ac.uk/pride/link/to/identification");
    sm.setSpectraRef("ms_run[2]:index=7|ms_run[2]:index=9");
    sm.setSearchEngine("[MS, MS:1001477, SpectraST,]");
    sm.setBestSearchEngineScore("[MS, MS:1001419, SpectraST:discriminant score F, 0.7]");
    sm.setModifications("CHEMMOD:+Na-H");

//    reference factory.addOptionalColumn(SmallMoleculeColumn.SEARCH_ENGINE_SCORE, metadata.getMsRunMap().get(1));
    sm.setSearchEngineScore(metadata.getMsRunMap().get(1), "[MS,MS:1001171,Mascot score,50]|[MS,MS:1001155,Sequest:xcorr,2]");

//  reference factory.addAbundanceOptionalColumn(metadata.getAssayMap().get(1));
    sm.setAbundanceColumn(metadata.getAssayMap().get(1), "12.3");


//  reference factory.addOptionalColumn(metadata.getMsRunMap().get(1), "my_value", String.class);
    sm.setOptionColumn(metadata.getMsRunMap().get(1), "my_value", "Tom");

    CVParam param = new CVParam("MS", "MS:1002217", "decoy peptide", null);
//  reference factory.addOptionalColumn(param, String.class);
    sm.setOptionColumn(param, "111");

    return sm;
}

public void testTabFile() throws Exception {
    Metadata metadata = getMetadata();
    MZTabColumnFactory factory = getSMH(metadata);
    SmallMolecule record = getRecord(metadata, factory);

    MZTabFile tabFile = new MZTabFile(metadata);
    tabFile.setSmallMoleculeColumnFactory(factory);
    tabFile.addSmallMolecule(record);
    
    BufferedOutputStream stream =new BufferedOutputStream(new FileOutputStream("E:\\lipidomics\\20100210\\test.mztab.txt"));
    tabFile.printMZTab(stream);
    stream.close();
}
*/
  private void createN15MassList(){

/*    try {
      Vector excelContent = QuantificationThread.parseQuantExcelFile("E:\\Dancy\\20131217\\PE Values.xls", 0f, 0f, 2, 1, true, 0f, 0f, 0f, 0f);
      LinkedHashMap<String,Integer> classSequence = (LinkedHashMap<String,Integer>)excelContent.get(0);
      Hashtable<String,Vector<String>> analyteSequence = (Hashtable<String,Vector<String>>)excelContent.get(1);
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)excelContent.get(2);
      String quantFile = "E:\\Dancy\\20131217\\PE Values_N15.xls";    
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(quantFile));
      HSSFWorkbook resultWorkbook = new HSSFWorkbook();
      HSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
      for (String className : classSequence.keySet()){
        HSSFSheet resultSheet = resultWorkbook.createSheet(className);
        int rowCountOut = 3;
        HSSFRow outRow = resultSheet.createRow(rowCountOut);
        HSSFCell label = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("Name");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("dbs");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("C");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("H");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("O");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("P");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("N");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(8,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("Nn");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(9,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("M");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("mass(form[-H] name[-H])");
        label.setCellStyle(headerStyle);
        label = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);
        label.setCellValue("tR (min)");
        label.setCellStyle(headerStyle);
        for (String analyteName : analyteSequence.get(className)){
          for (String mod:quantObjects.get(className).get(analyteName).keySet()){
            QuantVO quant = quantObjects.get(className).get(analyteName).get(mod);
            rowCountOut++;
            String formula = quant.getAnalyteFormula();
            StringTokenizer tokenizer = new StringTokenizer(formula," ");
            int amountC = 0;
            int amountH = 0;
            int amountO = 0;
            int amountP = 0;
            int amountN = 0;
            while (tokenizer.hasMoreTokens()){
              String elPlusNumber = tokenizer.nextToken().trim();
              String number = "";
              if (elPlusNumber.startsWith("C")||elPlusNumber.startsWith("H")||elPlusNumber.startsWith("O")||
                  elPlusNumber.startsWith("P")||elPlusNumber.startsWith("N")){
                number = elPlusNumber.substring(1);
              }
              if (number==null || number.length()==0) continue;
              try{
                int amount = Integer.parseInt(number);
                if (elPlusNumber.startsWith("C")) amountC = amount;
                else if (elPlusNumber.startsWith("H")) amountH = amount;
                else if (elPlusNumber.startsWith("O")) amountO = amount;
                else if (elPlusNumber.startsWith("P")) amountP = amount;
                else if (elPlusNumber.startsWith("N")) amountN = amount;
              }catch(NumberFormatException nfx){
                
              }
            }
            HSSFRow firstRow = resultSheet.createRow(rowCountOut);
            rowCountOut++;
            HSSFRow secondRow = resultSheet.createRow(rowCountOut);
            HSSFCell cellOut = firstRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
            cellOut.setCellValue(quant.getAnalyteName());
            cellOut = secondRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
            cellOut.setCellValue("N15"+quant.getAnalyteName());
            cellOut = firstRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
            cellOut.setCellValue(":");
            cellOut = secondRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
            cellOut.setCellValue(":");
            cellOut = firstRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(quant.getDbs());
            cellOut = secondRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(quant.getDbs());
            cellOut = firstRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountC);
            cellOut = secondRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountC);
            cellOut = firstRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountH);
            cellOut = secondRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountH);
            cellOut = firstRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountO);
            cellOut = secondRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountO);
            cellOut = firstRow.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountP);
            cellOut = secondRow.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountP);
            cellOut = firstRow.createCell(7,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountN);
            cellOut = secondRow.createCell(7,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(0);
            cellOut = firstRow.createCell(8,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(0);
            cellOut = secondRow.createCell(8,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(amountN);
            cellOut = firstRow.createCell(9,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(quant.getAnalyteMass()-1.007276);
            cellOut = secondRow.createCell(9,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(quant.getAnalyteMass()-1.007276+0.997);
            cellOut = firstRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(quant.getAnalyteMass());
            cellOut = secondRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
            cellOut.setCellValue(quant.getAnalyteMass()+0.997);
            
          }
        }
      }
      resultWorkbook.write(out);
      out.close();

    }
    catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch (SpectrummillParserException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch (ExcelInputFileException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch (ChemicalFormulaException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
   */ 
  }

  private void readMSnIdentification(){
    //String filePath = "C:\\Sphingolipids\\dataStandardsPrelim\\20190111\\positive\\01092018 CER and SPH Mixes positive-10uM MIX 1_Cer_test_pos.xlsx";
    String filePath = "C:\\data\\BiologicalExperiment\\Orbitrap_CID\\positive\\002_liver2-1_Orbitrap_CID_pos_positive.xlsx";
    try {
      Hashtable<String,Boolean> showMods = new  Hashtable<String,Boolean>();
      QuantificationResult result = LDAResultReader.readResultFile(filePath, showMods);
//      for (String key : result.getIdentifications().keySet()){
//        for (LipidParameterSet set : result.getIdentifications().get(key)){
//          if (!(set instanceof LipidomicsMSnSet)) continue;
//          LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
//          if (key.equalsIgnoreCase("DG") && set.getNameStringWithoutRt().equalsIgnoreCase("IS40:10")){
//            Hashtable<String,Hashtable<Integer,Integer>> posDef = msn.getPositionDefinition();
//            Hashtable<String,Hashtable<Integer,Vector<IntensityPositionVO>>> posEv = msn.getPositionEvidence();
//            for (String combi : posDef.keySet()){
//              System.out.println(combi);
//              Hashtable<Integer,Integer> pd = posDef.get(combi);
//              Hashtable<Integer,Vector<IntensityPositionVO>> pev = posEv.get(combi);
//              for (Integer order : pd.keySet()){
//                System.out.println(order+": "+pd.get(order)+";"+pev.get(pd.get(order)+1).size());
//                if (!pev.containsKey(pd.get(order)+1)) continue;
//                for (IntensityPositionVO intPos : pev.get(pd.get(order)+1)){
//                  System.out.println(intPos.getReadableRuleInterpretation());
//                }
//              }
//            }
//          }
//        }
//      }
    }
    catch (ExcelInputFileException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  private void calculateErrorPropagatedStdev(){
    double[] numerators = new double[3];
    numerators[0] = 127145.872309953d;
    numerators[1] = 148517.993889645d;
    numerators[2] = 164652.100566597d;
    double[] denominators = new double[3];
    denominators[0] = 131560.439614663d;
    denominators[1] = 126301.468674693d;
    denominators[2] = 137584.243909604d;
    double stdev = Calculator.calculateRatioStdevErrorPropagated(numerators,denominators,false);
    System.out.println("StDev: "+stdev);
  }
  
  private void sumVarianceErrorPropagated(){
    Vector<Double> values = new Vector<Double>();
    values.add(4.95572466632494d);
    values.add(4.11858097188603d);
    values.add(4.48943333538745d);
    values.add(5.84853127770684d);
    values.add(3.57847708039813d);
    values.add(3.78127478735857d);
    values.add(8.73964926013728d);
    System.out.println("Var: "+Calculator.calculateSumStdevErrorPropagated(values)/Math.sqrt(7));
  }
  
  
  /* Ergebnisse aus Levenberg-Marquadt
   * A = 0.084816925
   * B = -5.1918063
   * C = -0.06733154
   * D = -0.7511424
   * E = 102.1479
   */
  private float calculateRetentionTime(int cAtoms, int dbs){
    float a = 0.084816925f;
    float b = -5.1918063f;
    float c = -0.06733154f;
    float d = -0.7511424f;
    float e = 102.1479f;
    return (float)(Math.pow(cAtoms, 2)*a+cAtoms*b+Math.pow(dbs, 2)*c+dbs*d+e);
  }
  
  
  
  
  private float[] calculateRetentionTime2(int cAtoms, int dbs){
    float[] equResults = new  float[2];
//    double a = 16.045525d;
//    double b = 0.027311193d;
//    double c = 23.697464d;
//    double d = 0.07704608f;
//    double e = -39.4736f;
//    double a = -24073.234;
//    double b = -0.17241788;
//    double c = 7.282743;
//    double d = 0.18594968;
//    double e = 39.021255;
    /****    double a = 6.002452;
    double b = -0.10829978;
    double c = -6.1439137;
    double d = -0.38661823;
    double e = 18.997349;
    double a = 0.45586327;
    double b = 1.9720105;
    double c = 1.50897828E10;
    double d = 17.040405;*/
    
    /****double a = 0.68126184;
    double b = 0.052865483;
    double c = 16.029482;
    double d = 0.05046112;
    double e = 15.657554;*/
    
/***    double a = 0.11315229;
    double b = -6.2656374;
    double c = 0.16806002;
    double d = 6.828068;
    double e = -0.2551649;
    double f = 106.66271;
    
    double a = 0.98916054;
    double b = 0.073833525;
    double c = 12.258491;
    double d = 0.22200058;
    double e = 1.4969523;
    */
    
/*    double a = 0.5637313;
    double b = 7.674056;
    double c = 0.13044955;
    double d = 0.002735732;
    double e = 4.0643334;*/
    
/*    double a = 0.48954707;
    double b = 2.6126845;
    double c = 0.015216565;
    double d = 15.559856;*/
    
/*    double a = 21.013233;
    double b = 0.14042547;
    double c = 13.336806;
    double d = 0.05277403;
    double e = 0.006699581;
    double f = -17.918058;*/
    
/*     double a = 19.011599;
    double b = 0.025909346;
    double c = 9.3087225;
    double d = 0.061780863;
    double e = 0.03233426;
    double f = -1.2009382;
    
    double a = 25.59345;
    double b = 0.09342736;
    double c = 12.778386;
    double d = 0.06821823;
    double e = 0.0037436788;
    double f = -13.481577;*/
    
/*    double a = 21.324095;
    double b = 0.11970881;
    double c = 13.872888;
    double d = 0.054845177;
    double e = 0.006127314;
    double f = -15.319756;*/
    
    double a = 1.9967897;
    double b = 404.4479;
    double c = 14.31973;
    double d = 0.038756903;
    double e = 0.009066881;
    double f = 33.258083;
    

    
    //double equationResult = a*Math.exp(b*cAtoms)+c*Math.exp(-1d*d*dbs)+e;
    //double equationResult = a*cAtoms+b*Math.exp(d*cAtoms-c*dbs)+e;
    //double equationResult = a*cAtoms-b*(1f-c*cAtoms)*dbs+d;
    double equationResult = a*Math.log(b*cAtoms)+c*Math.exp(e*cAtoms-d*dbs)+f;
    //double equationResult = a*(1-b/cAtoms)+c*Math.exp(e*cAtoms-d*dbs)+f;
    
/*    double a2 = 1.3564086;
    double b2 = 274.8564;
    double c2 = 0.6016927;
    double d2 = 14.142465;
    double e2 = 0.050176375;
    double f2 = 0.0056969444;
    double g2 = 57.000843;
    
    double a2 = 0.12856463;
    double b2 = 4394.7827;
    double c2 = 0.7734002;
    double d2 = 14.221329;
    double e2 = 0.05237847;
    double f2 = 0.0056461473;
    double g2 = 50.38823;*/
    
    double a2 = 0.014360017;
    double b2 = 15423.295;
    double c2 = 0.45447597;
    double d2 = 12.373343;
    double e2 = 0.05145469;
    double f2 = 0.0075287553;
    double g2 = 44.85662;
    
    double equationResult2 = a2*(1-b2*Math.pow(cAtoms, -1f*c2))+d2*Math.exp(f2*cAtoms-e2*dbs)+g2;
    
    equResults[0] = (float)equationResult;
    equResults[1] = (float)equationResult2;
    
    return equResults;
  }

  
  
  private LevenbergMarquardtOptimizer leastSquareFittingFor2Variables() throws LMException{
    float[][] values = new float[15][2];
    values[0] = new float[]{40,6};
    values[1] = new float[]{40,5};
    values[2] = new float[]{40,5};
    values[3] = new float[]{40,4};
    values[4] = new float[]{38,6};
    values[5] = new float[]{38,6};
    values[6] = new float[]{38,5};
    values[7] = new float[]{38,4};
    values[8] = new float[]{38,3};
    values[9] = new float[]{38,3};
    values[10] = new float[]{36,4};
    values[11] = new float[]{36,3};
    values[12] = new float[]{36,2};
    values[13] = new float[]{34,2};
    values[14] = new float[]{34,1};
    
    
    float[] rtValues = new float[15];
    rtValues[0] = 23.35f;
    rtValues[1] = 24.01f;
    rtValues[2] = 24.85f;
    rtValues[3] = 25.42f;
    rtValues[4] = 20.78f;
    rtValues[5] = 20.23f;
    rtValues[6] = 21.85f;
    rtValues[7] = 23.87f;
    rtValues[8] = 24.55f;
    rtValues[9] = 25.44f;
    rtValues[10] = 21.28f;
    rtValues[11] = 21.91f;
    rtValues[12] = 23.71f;
    rtValues[13] = 21.10f;
    rtValues[14] = 22.94f;
    
    /*    float[][] values = new float[12][2];
    values[0] = new float[]{48,0};
    values[1] = new float[]{50,0};
    values[2] = new float[]{50,1};
    values[3] = new float[]{50,3};
    values[4] = new float[]{52,0};
    values[5] = new float[]{52,1};
    values[6] = new float[]{52,2};
    values[7] = new float[]{54,0};
    values[8] = new float[]{54,1};
    values[9] = new float[]{54,2};
    values[10] = new float[]{54,6};
    values[11] = new float[]{60,0};
//    values[12] = new float[]{62,1};
    
    float[] rtValues = new float[12];
    rtValues[0] = 40.85f;
    rtValues[1] = 41.82f;
    rtValues[2] = 39.96f;
    rtValues[3] = 39.21f;
    rtValues[4] = 42.77f;
    rtValues[5] = 41.95f;
    rtValues[6] = 41.07f;
    rtValues[7] = 43.65f;
    rtValues[8] = 42.82f;
    rtValues[9] = 42.03f;
    rtValues[10] = 39.10f;
    rtValues[11] = 46.32f;*/
//    rtValues[12][0] = 46.39f;
    
    LevenbergMarquardtOptimizer lmOptimizer = new LMQuadraticTwoVariables(values,rtValues,null);
    lmOptimizer.fit();
    
    System.out.println("Lambda "+lmOptimizer.getResultLambda());

//    residues = calculateResidues(resultVector,values,paramsVector);
//    for (int i=0;i!=residues.A.length;i++){
//      System.out.println(i+";"+residues.A[i][0]);
//    }

    System.out.println("Params:");
    FloatMatrix paramsVector = lmOptimizer.getResultParams();
    for (int i=0; i!=paramsVector.A.length; i++){
      System.out.println(paramsVector.A[i][0]);
    }
    System.out.println("Mean Deviation: "+lmOptimizer.getMeanDeviation());
    System.out.println("------------------------");
    return lmOptimizer;
  }
  
  
  private void calculateRetentionTime(){
//    try {

//      LevenbergMarquardtOptimizer lmOptimizer = leastSquareFittingFor2Variables();
      for (int cAtoms=26; cAtoms<70; cAtoms +=2){
        for (int dbs=0; dbs<15; dbs++){
          float[] results = calculateRetentionTime2(cAtoms,dbs);
          System.out.println(cAtoms+":"+dbs+": "+results[0]+" ; "+results[1]);
          float[] values = new float[2];
          values[0] = cAtoms;
          values[1] = dbs;
//          System.out.println(cAtoms+":"+dbs+": "+lmOptimizer.calculateFitValue(values));
        }
      }
//    }
//    catch (LMException e) {
//      // TODO Auto-generated catch block
//      e.printStackTrace();
//    }       
        
  }
  
  private void paintRtDependencies(){
/*    try {
      Hashtable<String,Vector<LipidParameterSet>> hits =  LipidDataAnalyzer.readResultFile("E:\\lipidomics\\20100210\\20100126_TAG-34_Massenliste-TAG_mit IS.xls", new Hashtable<String,Boolean>()).getIdentifications();
      String mod = "NH4";
      for (String className : hits.keySet()){
        Vector<LipidParameterSet> classHits = hits.get(className);
        String ruleName = "TG_"+mod;
        Pattern cAtomsPattern =  Pattern.compile(RulesContainer.getCAtomsFromNamePattern(ruleName));
        Pattern dbsPattern =  Pattern.compile(RulesContainer.getDoubleBondsFromNamePattern(ruleName));
        Range cAtomsRange = new Range(Integer.MAX_VALUE,0f);
        Hashtable<Integer,Range> dbsRanges = new Hashtable<Integer,Range>();
        // find the MSn identifications and group them according to their # C atoms and # double bonds
        Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> paramsOrdered = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
        for (LipidParameterSet set : classHits){
          String analyteName = set.getNameStringWithoutRt();
          Matcher cAtomsMatcher = cAtomsPattern.matcher(analyteName);
          if (!cAtomsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_CATOMS_PARSE+" pattern \""+RulesContainer.getCAtomsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
          int cAtoms = Integer.parseInt(cAtomsMatcher.group(1));
          Matcher dbsMatcher = dbsPattern.matcher(analyteName);
          if (!dbsMatcher.matches()) throw new RulesException("The analyte "+analyteName+" does not match the "+FragRuleParser.GENERAL_DBOND_PARSE+" pattern \""+RulesContainer.getDoubleBondsFromNamePattern(ruleName)+"\" of the class "+ruleName+"!");
          int dbs = Integer.parseInt(dbsMatcher.group(1));
          if (cAtoms<cAtomsRange.getStart()) cAtomsRange = new Range(cAtoms,cAtomsRange.getStop());
          if (cAtoms>cAtomsRange.getStop()) cAtomsRange = new Range(cAtomsRange.getStart(),cAtoms);
          Range dbsRange = new Range(Integer.MAX_VALUE,0f);
          if (dbsRanges.containsKey(cAtoms)) dbsRange = dbsRanges.get(cAtoms);
          if (dbs<dbsRange.getStart()) dbsRange = new Range(dbs,dbsRange.getStop());
          if (dbs>dbsRange.getStop()) dbsRange = new Range(dbsRange.getStart(),dbs);
          dbsRanges.put(cAtoms, dbsRange);
          
          String rt = set.getRt();
//          Hashtable<String,LipidParameterSet> hitsOfAnalyte = unprocessed.get(analyteName);
//          for (String rt : hitsOfAnalyte.keySet()){
//            LipidParameterSet params = hitsOfAnalyte.get(rt);
//          if ((analyteName.equalsIgnoreCase("50:1")&&rt.equalsIgnoreCase("39.96")) ||
//              (analyteName.equalsIgnoreCase("50:3")&&rt.equalsIgnoreCase("38.14")) ||
//              (analyteName.equalsIgnoreCase("52:2")&&rt.equalsIgnoreCase("40.14")) ||
//              (analyteName.equalsIgnoreCase("54:1")&&rt.equalsIgnoreCase("41.95")) ||
//              (analyteName.equalsIgnoreCase("54:2")&&rt.equalsIgnoreCase("41.10")) ||
//              (analyteName.equalsIgnoreCase("32:0")))
//            continue;
//            if (set instanceof LipidomicsMSnSet && !analyteName.startsWith("IS")){
              Hashtable<Integer,Hashtable<String,LipidParameterSet>> paramsSameC = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
              if (paramsOrdered.containsKey(cAtoms)) paramsSameC = paramsOrdered.get(cAtoms);
              Hashtable<String,LipidParameterSet> paramsSameDbs = new  Hashtable<String,LipidParameterSet>();
              if (paramsSameC.containsKey(dbs)) paramsSameDbs = paramsSameC.get(dbs);
              paramsSameDbs.put(rt, set);
              paramsSameC.put(dbs,paramsSameDbs);
              paramsOrdered.put(cAtoms, paramsSameC);
//            }
//          }
        }
        
        System.out.println("C-Range: "+Math.round(cAtomsRange.getStart())+";"+Math.round(cAtomsRange.getStop()));
//        int lowestDbs = Integer.MAX_VALUE;
//        int highestDbs = 0;
        for (Integer cAtoms : dbsRanges.keySet()){
          Range dbsRange = dbsRanges.get(cAtoms);
          System.out.println("Dbs-range "+cAtoms+": "+Math.round(dbsRange.getStart())+";"+Math.round(dbsRange.getStop()));
//          if (dbsRange.getStart()<lowestDbs) lowestDbs = Math.round(dbsRange.getStart());
//          if (dbsRange.getStop()>highestDbs) highestDbs = Math.round(dbsRange.getStop());
        }
        String dir = "E:\\lipidomicsMS2\\20140220\\rtDependency\\";
        Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>> dbsSorted = new Hashtable<Integer,Hashtable<Integer,Hashtable<String,LipidParameterSet>>>();
        
        String bothVarFile = dir+"bothVar.txt";
        String bothVar = "";
        
        for (Integer cAtoms : paramsOrdered.keySet()){
          String fileName = dir+"C"+String.valueOf(cAtoms)+".png";
          String fileNameText = dir+"C"+String.valueOf(cAtoms)+".txt";
          String text = "";
          XYSeries series = new XYSeries("double bonds");
          
          for (Integer dbs : paramsOrdered.get(cAtoms).keySet()){
            Hashtable<String,LipidParameterSet> sameRt = paramsOrdered.get(cAtoms).get(dbs);
            for (String rt : sameRt.keySet()){
              series.add(dbs, Double.valueOf(rt));
              text+=dbs+"\t"+rt+"\n";
              bothVar+=cAtoms+"\t"+dbs+"\t"+rt+"\n";
              System.out.println("MS2-Hit: "+cAtoms+":"+dbs+"_"+rt);
            }
            Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameDbs = new Hashtable<Integer,Hashtable<String,LipidParameterSet>>();
            if (dbsSorted.containsKey(dbs)) sameDbs = dbsSorted.get(dbs);
            sameDbs.put(cAtoms, sameRt);
            dbsSorted.put(dbs, sameDbs);
          }
          XYSeriesCollection coll = new XYSeriesCollection(series);
          JFreeChart chart = ChartFactory.createScatterPlot("C="+cAtoms, "dbs", "RT [min]", coll, PlotOrientation.VERTICAL,
            true, true, false);
//          ((Zoomable)chart.getPlot()).zoomRangeAxes(1.5, 2.0d, null, null);
          BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(fileName));
          ImageIO.write(chart.createBufferedImage(1000, 700), "PNG", stream);
          stream.close();
          stream = new BufferedOutputStream(new FileOutputStream(fileNameText));
          stream.write(text.getBytes());
          stream.close();
        }
        for (Integer dbs : dbsSorted.keySet()){
          String fileName = dir+"dbs"+String.valueOf(dbs)+".png";
          XYSeries series = new XYSeries("C Atoms");

          Hashtable<Integer,Hashtable<String,LipidParameterSet>> sameDbs = dbsSorted.get(dbs);
          if (sameDbs.size()<2) continue;
          for (Integer cAtoms : sameDbs.keySet()){
            Hashtable<String,LipidParameterSet> sameRt = sameDbs.get(cAtoms);
            for (String rt : sameRt.keySet()){
              series.add(cAtoms, Double.valueOf(rt));
            }
          }
          XYSeriesCollection coll = new XYSeriesCollection(series);
          JFreeChart chart = ChartFactory.createScatterPlot("dbs="+dbs, "C-Atoms", "RT [min]", coll, PlotOrientation.VERTICAL,
            true, true, false);
          BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(fileName));
          ImageIO.write(chart.createBufferedImage(1000, 700), "PNG", stream);
          stream.close();

        }
        BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(bothVarFile));
        stream.write(bothVar.getBytes());
        stream.close();
      }
    }
    catch (Exception e) {
      e.printStackTrace();
    }*/
  }
  
/*  private void paintPCAOfMostProminentComponents(){
    try {
      Hashtable<Integer,Hashtable<String,Double>> loadings = readPCALoadings();
      Hashtable<String,Hashtable<String,Double>> data = readPCAData();
      Hashtable<String,String> mainSpecies = new Hashtable<String,String>();
      
      XYSeries series1 = new XYSeries("WT-FED");
      XYSeries series2 = new XYSeries("KO-FED");
      XYSeries series3 = new XYSeries("WT-FAS");
      XYSeries series4 = new XYSeries("KO-FAS");*/
/****      mainSpecies.put("52:2", "52:2");
      mainSpecies.put("52:3", "52:3");
      mainSpecies.put("52:4", "52:4");
      mainSpecies.put("52:5", "52:5");
      mainSpecies.put("54:4", "54:4");
      mainSpecies.put("54:5", "54:5");
      mainSpecies.put("54:6", "54:6");
      mainSpecies.put("56:6", "56:6");
      mainSpecies.put("56:7", "56:7");
      mainSpecies.put("56:8", "56:8");*/
/*      
      mainSpecies.put("34:1", "34:1");
      mainSpecies.put("34:2", "34:2");
      mainSpecies.put("36:2", "36:2");
      mainSpecies.put("36:4", "36:4");
      mainSpecies.put("38:4", "38:4");
      mainSpecies.put("38:6", "38:6");
      
      Vector<Double> wtFedPC1 = new Vector<Double>();
      Vector<Double> wtFedPC2 = new Vector<Double>();
      Vector<Double> koFedPC1 = new Vector<Double>();
      Vector<Double> koFedPC2 = new Vector<Double>();
      Vector<Double> wtFasPC1 = new Vector<Double>();
      Vector<Double> wtFasPC2 = new Vector<Double>();
      Vector<Double> koFasPC1 = new Vector<Double>();
      Vector<Double> koFasPC2 = new Vector<Double>();
      
      double sumPC1 = 0d;
      double sumPC2 = 0d;
      for (String expName : data.keySet()){
        Hashtable<String,Double> perExp = data.get(expName);
        double pc1 = 0d;
        double pc2 = 0d;
        for (String analyte : perExp.keySet()){
          if (!mainSpecies.containsKey(analyte)) continue;
          pc1 += perExp.get(analyte)*loadings.get(1).get(analyte);
          pc2 += perExp.get(analyte)*loadings.get(2).get(analyte);
        }
        pc1 = pc1/10d;
        pc2 = pc2/10d;
        sumPC1 += pc1;
        sumPC2 += pc2;
        if (expName.equalsIgnoreCase("46")||expName.equalsIgnoreCase("47")||expName.equalsIgnoreCase("48")){
          wtFedPC1.add(pc1);
          wtFedPC2.add(pc2);
        } else if (expName.equalsIgnoreCase("49")||expName.equalsIgnoreCase("50")||expName.equalsIgnoreCase("51")){
          koFedPC1.add(pc1);
          koFedPC2.add(pc2);
        } else if (expName.equalsIgnoreCase("52")||expName.equalsIgnoreCase("53")||expName.equalsIgnoreCase("54")){
          wtFasPC1.add(pc1);
          wtFasPC2.add(pc2);
        } else if (expName.equalsIgnoreCase("55")||expName.equalsIgnoreCase("56")||expName.equalsIgnoreCase("57")){
          koFasPC1.add(pc1);
          koFasPC2.add(pc2);
        }
      }
      String pc1String = "c(";
      String pc2String = "c(";
      for (int i=0; i!=wtFedPC1.size();i++){
        if (i!=0){
          pc1String += ",";
          pc2String += ",";
        }
        pc1String += String.valueOf(wtFedPC1.get(i)-sumPC1/((double)data.size()));
        pc2String += String.valueOf(wtFedPC2.get(i)-sumPC2/((double)data.size()));
      }
      for (int i=0; i!=koFedPC1.size();i++){
        pc1String += ","+String.valueOf(koFedPC1.get(i)-sumPC1/((double)data.size()));
        pc2String += ","+String.valueOf(koFedPC2.get(i)-sumPC2/((double)data.size()));
      }
      for (int i=0; i!=wtFasPC1.size();i++){
        pc1String += ","+String.valueOf(wtFasPC1.get(i)-sumPC1/((double)data.size()));
        pc2String += ","+String.valueOf(wtFasPC2.get(i)-sumPC2/((double)data.size()));
      }
      for (int i=0; i!=koFasPC1.size();i++){
        pc1String += ","+String.valueOf(koFasPC1.get(i)-sumPC1/((double)data.size()));
        pc2String += ","+String.valueOf(koFasPC2.get(i)-sumPC2/((double)data.size()));
      }
      pc1String += ")";
      pc2String += ")";
      System.out.println(pc1String);
      System.out.println(pc2String);
      
      for (int i=0;i!=wtFedPC1.size();i++){
        series1.add(wtFedPC1.get(i)-sumPC1/((double)data.size()),wtFedPC2.get(i)-sumPC2/((double)data.size()));
        series2.add(koFedPC1.get(i)-sumPC1/((double)data.size()),koFedPC2.get(i)-sumPC2/((double)data.size()));
        series3.add(wtFasPC1.get(i)-sumPC1/((double)data.size()),wtFasPC2.get(i)-sumPC2/((double)data.size()));
        series4.add(koFasPC1.get(i)-sumPC1/((double)data.size()),koFasPC2.get(i)-sumPC2/((double)data.size()));
      }
      
      XYSeriesCollection coll = new XYSeriesCollection(series1);
      coll.addSeries(series2);
      coll.addSeries(series3);
      coll.addSeries(series4);
      JFreeChart chart = ChartFactory.createScatterPlot("ATGL TG PCA", "PC1", "PC2", coll, PlotOrientation.VERTICAL,
        true, true, false);
//      ((Zoomable)chart.getPlot()).zoomRangeAxes(1.5, 2.0d, null, null);*/
/***      BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream("L:\\PaperNewJournal\\PCA_TG.png"));
      ImageIO.write(chart.createBufferedImage(1000, 700), "PNG", stream);
      stream.close();*/
/*    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
  }*/
  
  private void findPCAMostContributingComponents(){
    try {
      Hashtable<Integer,Hashtable<String,Double>> loadings = readPCALoadings();
      Hashtable<String,Hashtable<String,Double>> data = readPCAData();
      Vector<String> group1 = new Vector<String>();
      group1.add("46");
      group1.add("47");
      group1.add("48");
      Vector<String> group2 = new Vector<String>();
      group2.add("49");
      group2.add("50");
      group2.add("51");
      
      System.out.println("------------------- PC1 -------------------");
      Hashtable<String,Double> loads = loadings.get(1);
      calculatePCADifferences(loads,data,group1,group2);
      
      System.out.println("------------------- PC2 -------------------");
      loads = loadings.get(2);
      calculatePCADifferences(loads,data,group1,group2);
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  private void calculatePCADifferences(Hashtable<String,Double> loads,Hashtable<String,Hashtable<String,Double>> data,Vector<String> group1, Vector<String> group2){
    List<KeyValueVO> vos = new ArrayList<KeyValueVO>();
    for (String analyte : loads.keySet()){
      Vector<Double> values1 = new Vector<Double>();
      for (String exp : group1){
        if (data.containsKey(exp) && data.get(exp).containsKey(analyte)) values1.add(data.get(exp).get(analyte));
      }
      double mean1 = Calculator.mean(values1);
      Vector<Double> values2 = new Vector<Double>();
      for (String exp : group2){
        if (data.containsKey(exp) && data.get(exp).containsKey(analyte)) values2.add(data.get(exp).get(analyte));
      }
      double mean2 = Calculator.mean(values2);
      KeyValueVO vo = new KeyValueVO(analyte,(mean1-mean2)*loads.get(analyte));
      vos.add(vo);
    }
    Collections.sort(vos, new GeneralComparator("at.tugraz.genome.TestClass$KeyValueVO", "getValue", "java.lang.Double"));
    for (int i=(vos.size()-1);i!=-1;i--){
//    for (int i=0;i!=vos.size();i++){
      KeyValueVO vo = vos.get(i);
      System.out.println(vo.getName()+"; "+vo.getValue());
    }
  }
  
  private Hashtable<Integer,Hashtable<String,Double>> readPCALoadings() throws Exception{
    Hashtable<Integer,Hashtable<String,Double>> loadings = new Hashtable<Integer,Hashtable<String,Double>>();
    
//    LineNumberReader reader = new LineNumberReader(new FileReader("L:\\Paper ATGL\\PCA-Analysis\\PC_loadings.txt"));
    LineNumberReader reader = new LineNumberReader(new FileReader("L:\\PaperNewJournal\\secondSubmission\\PCA\\Lipidome_loadings.txt"));
    String line;
    while ((line = reader.readLine())!=null){
      StringTokenizer tokenizer = new StringTokenizer(line,"\t");
      if (tokenizer.countTokens()!=4)continue;
      String keyPCA1 = tokenizer.nextToken();//.substring(2);
      Double valuePCA1 = new Double(tokenizer.nextToken());
      String keyPCA2 = tokenizer.nextToken();//.substring(2);
      Double valuePCA2 = new Double(tokenizer.nextToken());
      Hashtable<String,Double> kvP1 = new Hashtable<String,Double>();
      Hashtable<String,Double> kvP2 = new Hashtable<String,Double>();
      if (loadings.containsKey(1)) kvP1 = loadings.get(1);
      if (loadings.containsKey(2)) kvP2 = loadings.get(2);
      kvP1.put(keyPCA1, valuePCA1);
      kvP2.put(keyPCA2, valuePCA2);
      loadings.put(1, kvP1);
      loadings.put(2, kvP2);
    }
    return loadings;
  }
  
  private Hashtable<String,Hashtable<String,Double>> readPCAData() throws Exception{
    Hashtable<String,Hashtable<String,Double>> data = new Hashtable<String,Hashtable<String,Double>>();
    LineNumberReader reader = new LineNumberReader(new FileReader("D:\\Development\\LipidomicNet\\ATGL\\lipidome.txt"));
    String line;
    int lineCount = 0;
    Hashtable<Integer,String> lookup = new Hashtable<Integer,String>();
    while ((line = reader.readLine())!=null){
      lineCount++;
      StringTokenizer tokenizer = new StringTokenizer(line,"\t");
      if (lineCount<2) continue;
      else if (lineCount==2){
        tokenizer.nextToken();
        int columnCount = 0;
        while(tokenizer.hasMoreTokens()){
          String expName = tokenizer.nextToken();
          lookup.put(columnCount, expName);
          columnCount++;
        }
      }else{
        String analyte = tokenizer.nextToken();
        tokenizer.nextToken();
        int columnCount = 0;
        while(tokenizer.hasMoreTokens()){
          Double value = new Double(tokenizer.nextToken());
          String columnName = lookup.get(columnCount);
          Hashtable<String,Double> values = new Hashtable<String,Double>();
          if (data.containsKey(columnName)) values = data.get(columnName);
          values.put(analyte, value);
          data.put(columnName, values);
          columnCount++;
        }        
      }
    }
    return data;
  }
  
  public class KeyValueVO{
    
    private String name_;
    private double value_;
    
    public KeyValueVO(String name, double value){
      this.name_ = name;
      this.value_ = value;
    }

    public String getName()
    {
      return name_;
    }

    public double getValue()
    {
      return value_;
    }
    
    
  }
    
  private void makeCommaSeparatedFiles(){
    try{
      String dir = "L:\\PaperNewJournal\\reducedPCA\\PC";
      File directory = new File(dir);
      for (File file : directory.listFiles()){
        if (file.isDirectory() || !file.getAbsolutePath().endsWith(".txt")) continue;
        LineNumberReader reader = new LineNumberReader(new FileReader(file));
        BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(file.getAbsolutePath().substring(0,file.getAbsolutePath().length()-3)+"tsv"));
        String line;
        while ((line = reader.readLine())!=null){
          StringTokenizer tokenizer = new StringTokenizer(line,";");
          String correctedLine = "";
          int count = 0;
          while (tokenizer.hasMoreTokens()){
            if (count!=0) correctedLine += "\t";
            correctedLine += tokenizer.nextToken().trim();
            count++;
          }
          correctedLine += "\n";
          out.write(correctedLine.getBytes());
        }
        out.close();
        reader.close();
      }
    } catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void doPostQuantProcessing(){
    try {
      Hashtable<String,Vector<LipidParameterSet>> results = LDAResultReader.readResultFile("E:\\lipidomicsMS2\\20140218\\neg\\030_Liver1_Neg_PI.xlsx", new Hashtable<String,Boolean>()).getIdentifications();
      //Hashtable<String,Vector<LipidParameterSet>> results = LipidDataAnalyzer.readResultFile("E:\\lipidomicsMS2\\20140220\\001_LD1_TG_all.xlsx", new Hashtable<String,Boolean>()).getIdentifications();
      //Hashtable<String,Vector<LipidParameterSet>> results = LipidDataAnalyzer.readResultFile("E:\\lipidomicsMS2\\testTGold\\20100126_TAG-34_TG_all.xlsx", new Hashtable<String,Boolean>()).getIdentifications();
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> ms2Removed = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
      
      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> forProcessing = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
      for (String className : results.keySet()){
        Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultsClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
        Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> ms2RemovedClass = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
        for (LipidParameterSet set : results.get(className)){
          String analyteName = set.getNameStringWithoutRt();
          String modName = set.getModificationName();
          String rt = set.getRt();
//          if (modName.equalsIgnoreCase("NH4") && (
//              (analyteName.equalsIgnoreCase("28:0")&&rt.equalsIgnoreCase("21.96"))||(analyteName.equalsIgnoreCase("34:0")&&rt.equalsIgnoreCase("29.30"))||
//              (analyteName.equalsIgnoreCase("36:0")&&rt.equalsIgnoreCase("31.33"))||(analyteName.equalsIgnoreCase("38:1")&&rt.equalsIgnoreCase("31.92"))||
//              (analyteName.equalsIgnoreCase("40:0")&&rt.equalsIgnoreCase("34.64"))||(analyteName.equalsIgnoreCase("40:1")&&rt.equalsIgnoreCase("33.64"))||
//              (analyteName.equalsIgnoreCase("42:0")&&rt.equalsIgnoreCase("36.15"))||(analyteName.equalsIgnoreCase("42:1")&&rt.equalsIgnoreCase("35.03"))||
//              (analyteName.equalsIgnoreCase("42:2")&&rt.equalsIgnoreCase("33.96"))||(analyteName.equalsIgnoreCase("44:0")&&rt.equalsIgnoreCase("37.50"))||
//              (analyteName.equalsIgnoreCase("44:1")&&rt.equalsIgnoreCase("36.45"))||(analyteName.equalsIgnoreCase("44:2")&&rt.equalsIgnoreCase("35.35"))||
//              (analyteName.equalsIgnoreCase("44:3")&&rt.equalsIgnoreCase("34.15"))||(analyteName.equalsIgnoreCase("46:0")&&rt.equalsIgnoreCase("38.72"))||
//              (analyteName.equalsIgnoreCase("46:1")&&rt.equalsIgnoreCase("37.75"))||(analyteName.equalsIgnoreCase("46:2")&&rt.equalsIgnoreCase("36.67"))||
//              (analyteName.equalsIgnoreCase("48:2")&&rt.equalsIgnoreCase("38.00"))||(analyteName.equalsIgnoreCase("48:3")&&rt.equalsIgnoreCase("36.92"))||
//              (analyteName.equalsIgnoreCase("50:2")&&rt.equalsIgnoreCase("39.04"))||(analyteName.equalsIgnoreCase("52:3")&&rt.equalsIgnoreCase("41.47"))||
//              (analyteName.equalsIgnoreCase("52:4")&&rt.equalsIgnoreCase("38.32"))||(analyteName.equalsIgnoreCase("52:5")&&rt.equalsIgnoreCase("37.75"))||
//              (analyteName.equalsIgnoreCase("54:3")&&rt.equalsIgnoreCase("40.27"))||(analyteName.equalsIgnoreCase("54:4")&&rt.equalsIgnoreCase("39.35"))||
//              (analyteName.equalsIgnoreCase("54:5")&&rt.equalsIgnoreCase("38.52"))||(analyteName.equalsIgnoreCase("56:0")&&rt.equalsIgnoreCase("43.65"))||
//              (analyteName.equalsIgnoreCase("56:1")&&rt.equalsIgnoreCase("42.85"))||(analyteName.equalsIgnoreCase("56:2")&&rt.equalsIgnoreCase("42.07"))||
//              (analyteName.equalsIgnoreCase("56:3")&&rt.equalsIgnoreCase("41.25"))||(analyteName.equalsIgnoreCase("56:4")&&rt.equalsIgnoreCase("40.43"))||
//              (analyteName.equalsIgnoreCase("56:5")&&rt.equalsIgnoreCase("39.59"))||(analyteName.equalsIgnoreCase("56:6")&&rt.equalsIgnoreCase("38.70"))||
//              (analyteName.equalsIgnoreCase("56:6")&&rt.equalsIgnoreCase("39.07"))||(analyteName.equalsIgnoreCase("58:0")&&rt.equalsIgnoreCase("44.60"))||
//              (analyteName.equalsIgnoreCase("58:1")&&rt.equalsIgnoreCase("43.77"))||(analyteName.equalsIgnoreCase("58:2")&&rt.equalsIgnoreCase("42.99"))||
//              (analyteName.equalsIgnoreCase("58:5")&&rt.equalsIgnoreCase("40.65"))||(analyteName.equalsIgnoreCase("58:6")&&rt.equalsIgnoreCase("39.85"))||
//              (analyteName.equalsIgnoreCase("60:1")&&rt.equalsIgnoreCase("44.67"))||(analyteName.equalsIgnoreCase("60:2")&&rt.equalsIgnoreCase("43.90"))||
//              (analyteName.equalsIgnoreCase("60:5")&&rt.equalsIgnoreCase("41.60"))||(analyteName.equalsIgnoreCase("62:0")&&rt.equalsIgnoreCase("46.35"))||
//              (analyteName.equalsIgnoreCase("62:2")&&rt.equalsIgnoreCase("44.71"))||(analyteName.equalsIgnoreCase("62:3")&&rt.equalsIgnoreCase("43.90"))||
//              (analyteName.equalsIgnoreCase("64:0")&&rt.equalsIgnoreCase("47.20"))||(analyteName.equalsIgnoreCase("64:1")&&rt.equalsIgnoreCase("46.39"))
//              (analyteName.equalsIgnoreCase("44:1")&&rt.equalsIgnoreCase("20.71"))||(analyteName.equalsIgnoreCase("44:2")&&rt.equalsIgnoreCase("19.75"))||
//              (analyteName.equalsIgnoreCase("44:3")&&rt.equalsIgnoreCase("19.82"))||(analyteName.equalsIgnoreCase("46:5")&&rt.equalsIgnoreCase("19.41"))||
//              (analyteName.equalsIgnoreCase("46:6")&&rt.equalsIgnoreCase("18.92"))||(analyteName.equalsIgnoreCase("48:4")&&rt.equalsIgnoreCase("20.39"))||
//              (analyteName.equalsIgnoreCase("48:6")&&rt.equalsIgnoreCase("19.89"))||(analyteName.equalsIgnoreCase("48:6")&&rt.equalsIgnoreCase("18.92"))||
//              (analyteName.equalsIgnoreCase("48:6")&&rt.equalsIgnoreCase("19.85"))||
//              (analyteName.equalsIgnoreCase("58:7")&&rt.equalsIgnoreCase("23.21"))||(analyteName.equalsIgnoreCase("58:8")&&rt.equalsIgnoreCase("22.53"))||
//              (analyteName.equalsIgnoreCase("58:10")&&rt.equalsIgnoreCase("21.53"))||(analyteName.equalsIgnoreCase("IS51:1")&&rt.equalsIgnoreCase("24.89"))
//              )){
//
//            putToHash(ms2RemovedClass,set,analyteName,modName,rt);
//          }else{
          
            putToHash(resultsClass,set,analyteName,modName,rt);
//          }
        }
        forProcessing.put(className, resultsClass);
        ms2Removed.put(className, ms2RemovedClass);
      }
      PostQuantificationProcessor processor = new PostQuantificationProcessor(forProcessing,ms2Removed,null);
      processor.processData();
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  private void putToHash(Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultsClass,LipidParameterSet set,String analyteName,String modName,String rt){
    Hashtable<String,Hashtable<String,LipidParameterSet>> resultsAnalyte = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
    if (resultsClass.containsKey(analyteName)) resultsAnalyte = resultsClass.get(analyteName);
    Hashtable<String,LipidParameterSet> resultsMod = new Hashtable<String,LipidParameterSet>();
    if (resultsAnalyte.containsKey(modName)) resultsMod = resultsAnalyte.get(modName);
    resultsMod.put(rt, set);
    resultsAnalyte.put(modName,resultsMod);
    resultsClass.put(analyteName,resultsAnalyte);

  }
  
  private void mergeLipidClassesAndCalcPercentage(){
    try{
      String directory = "L:\\PaperNewJournal\\data\\totalLipidContent";
      File dir = new File(directory);
      File[] toParse = dir.listFiles();
      Hashtable<String,Hashtable<String,Double>> values = new Hashtable<String,Hashtable<String,Double>>();
      Vector<String> analytesInSequence = new Vector<String>();
      Vector<String> expsInSequence = new Vector<String>();
      for (int i=0; i!=toParse.length; i++){
        File lipidClassFile = toParse[i];
        if (!lipidClassFile.getName().endsWith(".txt") || lipidClassFile.getName().equalsIgnoreCase("lipidome.txt")) continue;
        Vector<String> analytesOfFile =  parseResultsAndAddToExperimentHash(lipidClassFile,values,expsInSequence);
        analytesInSequence.addAll(analytesOfFile);
      }
      Hashtable<String,Hashtable<String,Double>> relativeValues = getRelativeToTotalSum(values);
      writeRelativeResultsToTxt(directory+"\\lipidome.txt",expsInSequence,analytesInSequence,relativeValues);
    }catch (Exception ex){
      
    }
  }
  
  private Vector<String> parseResultsAndAddToExperimentHash(File lipidClassFile, Hashtable<String,Hashtable<String,Double>> values, Vector<String> expsInSequence) throws Exception{
    String line;
    Vector<String> analytesInSequence = new Vector<String>();
    LineNumberReader reader = new LineNumberReader(new FileReader(lipidClassFile));
    int lineCount = 0;
    String lipidClass = "";
    double multiplicator = 1d;
    Hashtable<Integer,String> columnToExperiment = new Hashtable<Integer,String>();
    while ((line = reader.readLine()) != null) {
      if (lineCount==0) lipidClass = line.trim();
      else if (lineCount==1){
        String[] columns = line.split("\t");
        String unit = columns[0].substring(columns[0].indexOf("[")+1,columns[0].indexOf("]"));
        if (unit.equalsIgnoreCase("pmol")){
        } else if (unit.equalsIgnoreCase("nmol")){
          multiplicator = 1000d;
        }else{
          throw new Exception("The unit "+unit+" is not supported!");
        }
        for (int i=2; i<columns.length; i++){
          columnToExperiment.put(i, columns[i].trim());
        }
      }else{
        String[] columns = line.split("\t");
        if (columns.length!=(columnToExperiment.size()+2))continue;
        String speciesName = lipidClass+columns[0].trim();
        analytesInSequence.add(speciesName);
        for (int column=2; column<columns.length; column++){
          double value = Double.parseDouble(columns[column].trim())*multiplicator;
          String expName = columnToExperiment.get(column);
          Hashtable<String,Double> valuesOfExp = new Hashtable<String,Double>();
          if (values.containsKey(expName)) valuesOfExp = values.get(expName);
          else expsInSequence.add(expName);
          valuesOfExp.put(speciesName, value);
          values.put(expName, valuesOfExp);
        }
      }
      lineCount++;
    }
    return analytesInSequence;
  }
  
  private Hashtable<String,Hashtable<String,Double>> getRelativeToTotalSum(Hashtable<String,Hashtable<String,Double>> absolute){
    Hashtable<String,Hashtable<String,Double>> relative = new Hashtable<String,Hashtable<String,Double>>();
    for (String exp: absolute.keySet()){
      Hashtable<String,Double> analytes = absolute.get(exp);
      double sum = 0d;
      for (String anal : analytes.keySet()) sum += analytes.get(anal);
      Hashtable<String,Double> analytesRel = new Hashtable<String,Double>();
      for (String anal : analytes.keySet()){
        analytesRel.put(anal, (analytes.get(anal)*1000d)/sum);
      }
      relative.put(exp, analytesRel);
    }
    return relative;
  }
  
  private void writeRelativeResultsToTxt(String fileName, Vector<String> expsInSequence, Vector<String> analytesInSequence, Hashtable<String,Hashtable<String,Double>> values) throws Exception{
    String content = "";
    content += "lipidome\n";
    content += "total amount [\u2030]\t";
    for (String exp : expsInSequence) content += "\t"+exp;
    content += "\n";
    for (String anal : analytesInSequence){
      content += anal+"\tValue";
      for (String exp : expsInSequence){
        Hashtable<String,Double> valuesOfExp =  values.get(exp);
        if (!valuesOfExp.containsKey(anal)) continue;
        content += "\t"+Calculator.FormatNumberToString(valuesOfExp.get(anal), 4d);
      }
      content+="\n";
    }
    
    BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(fileName));
    out.write(content.getBytes());
    out.close();
  }
  
  private void generateAlkylLinkedMassList(){
    try{
      ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
      parser.parse();
      String quantFile = "E:\\lipidomicsMS2\\20140528\\O-PE_new.xlsx";
      
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(quantFile));
      XSSFWorkbook resultWorkbook = new XSSFWorkbook();
      XSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
      Sheet resultSheet = resultWorkbook.createSheet("O-PE");
      int rowCount = 4;
      Row outRow = resultSheet.createRow(rowCount);
      Cell label = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("Name");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("dbs");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("C");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("H");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("O");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("P");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("N");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(8,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("D");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(9,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("M");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("mass(form[-H] name[-H])");
//      label.setCellValue("mass(form[-H] name[-H])");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("tR (min)");
      label.setCellStyle(headerStyle);
      rowCount++;
      
      Hashtable<Integer,Integer> cDbsCombi = new Hashtable<Integer,Integer>();
      cDbsCombi.put(24, 0);
      cDbsCombi.put(25, 0);
      cDbsCombi.put(26, 2);
      cDbsCombi.put(27, 2);
      cDbsCombi.put(28, 4);
      cDbsCombi.put(29, 4);
      cDbsCombi.put(30, 4);
      cDbsCombi.put(31, 4);
      cDbsCombi.put(32, 4);
      cDbsCombi.put(33, 4);
      cDbsCombi.put(34, 6);
      cDbsCombi.put(35, 6);
      cDbsCombi.put(36, 8);
      cDbsCombi.put(37, 8);
      cDbsCombi.put(38, 8);
      cDbsCombi.put(39, 8);
      cDbsCombi.put(40, 8);
      cDbsCombi.put(41, 8);
      cDbsCombi.put(42, 8);
      cDbsCombi.put(43, 8);
      cDbsCombi.put(44, 8);
      cDbsCombi.put(45, 8);
      cDbsCombi.put(46, 10);
      cDbsCombi.put(47, 10);
      cDbsCombi.put(48, 10);
      
      for (int cAtoms=24; cAtoms<49; cAtoms+=1){
        int dbs = 0;
        int dbsMax = cDbsCombi.get(cAtoms)+1;
        while (dbs<dbsMax){
          int totalC = cAtoms+5;
//          int totalC = cAtoms+8;
          int hAtoms = totalC*2-2*dbs+2;
//          int hAtoms = totalC*2-2*dbs+4;
//          int hAtoms = totalC*2-2*dbs;
//          int hAtoms = totalC*2-2*dbs;
          int oAtoms = 7;
          int pAtoms = 1;
          int nAtoms = 1;
          int dAtoms = 0;
          double massNeutral = parser.getElementDetails("C").getMonoMass()*((double)totalC)+parser.getElementDetails("H").getMonoMass()*((double)hAtoms)+
              parser.getElementDetails("O").getMonoMass()*((double)oAtoms)+parser.getElementDetails("P").getMonoMass()*((double)pAtoms)+
              parser.getElementDetails("N").getMonoMass()*((double)nAtoms)+parser.getElementDetails("D").getMonoMass()*((double)dAtoms);
//          double massCharged = massNeutral+parser.getElementDetails("H").getMonoMass()+parser.getElementDetails("C").getMonoMass()+2d*parser.getElementDetails("O").getMonoMass();
          double massCharged = massNeutral-parser.getElementDetails("H").getMonoMass();
          
          outRow = resultSheet.createRow(rowCount);
          label = outRow.createCell(0,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(cAtoms);
          label = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
          label.setCellValue(":");       
          label = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(dbs);
          label = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(totalC);
          label = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(hAtoms);
          label = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(oAtoms);
          label = outRow.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(pAtoms);
          label = outRow.createCell(7,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(nAtoms);
          label = outRow.createCell(8,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(dAtoms);
          label = outRow.createCell(9,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(massNeutral);
          label = outRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(massCharged);
          dbs++;
          rowCount++;
        }
      }
      
      resultWorkbook.write(out);
      out.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }

  }
  
  private void quantifyPeak(){
    String[] chromPaths = StringUtils.getChromFilePaths("D:\\lipidomics\\20150626\\28.chrom");   
    System.out.println(chromPaths[1]);
    LipidomicsAnalyzer lAnalyzer;
    ElementConfigParser aaParser = Settings.getElementParser();
    try {
      lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
      setStandardParameters(lAnalyzer);
      double mz = 758.569432004812;
      Vector<Double> mustMatchProbabs = aaParser.calculateChemicalFormulaIntensityDistribution("C42 H80 O8 P1 N1", 2, false).get(0);
      Vector<Double> probabs = aaParser.calculateChemicalFormulaIntensityDistribution("C42 H80 O8 P1 N1", 3, false).get(0);

      Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> probes = lAnalyzer.processByMzProbabsAndPossibleRetentionTime((float)mz, 1, null, 0f, 0f, LipidomicsDefines.MINUTES, mustMatchProbabs, probabs, 1, false);
      
      for (Integer key1 : probes.keySet()){
        Hashtable<Integer,Vector<CgProbe>> set1 = probes.get(key1);
        for (Integer key2 : set1.keySet()){
          Vector<CgProbe> set2 = set1.get(key2);
          for (CgProbe probe : set2) System.out.println(key1+";"+key2+";"+probe.Peak/60f+"  "+probe.Area);
        }
      }
      
      
//      for (Float rt : rts){
//        System.out.println("RT: "+rt);
//      }
      //MSnDebugVO debugInfo = analyzer.getDebugInfo();

    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

  }

  
  private void searchMSnSpectra(){
    String[] chromPaths = StringUtils.getChromFilePaths("E:\\Robin\\20140819\\Control\\10uL_21Mar14_L4440_Negative_MS2_32-97-32_300uLmin_01.chrom");   
    System.out.println(chromPaths[1]);
    LipidomicsAnalyzer lAnalyzer;
    ElementConfigParser aaParser = Settings.getElementParser();
    try {
      lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
      setStandardParameters(lAnalyzer);
      double mz = 638.40380859375;
      double tol = 0.015;
      long time = System.currentTimeMillis();
      MSnAnalyzer analyzer = new MSnAnalyzer("PC","HCOO",mz, tol,"22",0,0,"C30 H60 O8 P1 N1","HCOO",1,lAnalyzer,true,false);
      Vector<Float> rts = analyzer.getFoundMatchingSpectraTimes();
      Vector<Double> mustMatchProbabs = aaParser.calculateChemicalFormulaIntensityDistribution("C31 H61 O10 P1 N1", 2, false).get(0);
      Vector<Double> probabs = aaParser.calculateChemicalFormulaIntensityDistribution("C31 H61 O10 P1 N1", 3, false).get(0);
      long lasted = System.currentTimeMillis()-time;
      System.out.println("Lasted: "+lasted/1000l+" "+lasted%1000l);
      time = System.currentTimeMillis();
      Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> probes = lAnalyzer.processByMzProbabsAndPossibleRetentionTime((float)mz,1,rts, 0f, 0f, LipidomicsDefines.MINUTES, 
          mustMatchProbabs,probabs,1,false);
      lasted = System.currentTimeMillis()-time;
      System.out.println("Lasted: "+lasted/1000l+" "+lasted%1000l);
      
      for (Integer key1 : probes.keySet()){
        Hashtable<Integer,Vector<CgProbe>> set1 = probes.get(key1);
        for (Integer key2 : set1.keySet()){
          Vector<CgProbe> set2 = set1.get(key2);
          for (CgProbe probe : set2) System.out.println(key1+";"+key2+";"+probe.Peak/60f);
        }
      }
      
      
//      for (Float rt : rts){
//        System.out.println("RT: "+rt);
//      }
      //MSnDebugVO debugInfo = analyzer.getDebugInfo();

    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

  }
  
  private void mergeMzXMLFiles(){
    //MzXMLMergerForWaters merger = new MzXMLMergerForWaters("D:\\lipidomcsMS2\\20141112\\20141110_GNR_Triebl_Survey_40VCE_Test_01.mzXML",5);
    MzXMLMergerForWaters merger = new MzXMLMergerForWaters("D:\\lipidomics\\20150817\\20150731_CenPK0715_C13_0.5h.mzXML",2);
    try {
      merger.merge();
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  private void convertWiffFile(){
    String[] params  = new String[5];
    params[0] = "D:\\Development\\LipidDataAnalyzer\\installer\\windows\\x64\\msconvert.exe";
    params[1] = "--mzXML";
    params[2] = "L:\\Controlled_Exp_LDA2\\WiffTest\\Data20150105_Ex1_neg.wiff";
    params[3] = "-o";
    params[4] = "L:\\Controlled_Exp_LDA2\\WiffTest";
    for (String param : params){
      System.out.println(param);
    }
    Process process;
    try {
      process = Runtime.getRuntime().exec(params);
      process.waitFor();
      System.out.println("Finishedddddddddddddddddd");
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  private void evaluateExperiment3(){
    String baseDir = "D:\\Experiment3\\Orbitrap_CID\\positive\\";
//    String baseDir = "D:\\Experiment3\\Orbitrap_HCD\\positive\\";
//    String baseDir = "D:\\lipidomcsMS2\\20150527_Exp3\\CID\\";
//    String baseDir = "D:\\Experiment3\\QTOF\\positive_new-dontUse\\";
//    String baseDir = "D:\\Experiment3\\QTOF\\negative\\";
//    String baseDir = "D:\\lipidomcsMS2\\QTOF\\20151125_Exp3\\";
    
//    String baseDir = "D:\\Experiment3\\QTRAP\\positive\\";
//    String baseDir = "D:\\Experiment3\\Orbitrap_HCD\\negative\\";
    try{
      //BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(baseDir+"Exp3-QTRAP-result_pos.xlsx"));
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(baseDir+"Exp3-CID-result_pos.xlsx"));
      Workbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      CellStyle rightStyle = getRightStyle(resultWorkbook);
      Sheet sheet = resultWorkbook.createSheet("Result");
      int rowCount = 0;
      Row row = sheet.createRow(rowCount);
      rowCount++;
      Row secondRow = sheet.createRow(rowCount);
      rowCount++;
      int fileColumn = 0;
      setHeaderCells(sheet,row,secondRow,headerStyle,"File","",fileColumn);
      //positive Columns
      int pgColumn = 2;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PG","18:1/18:1-18:0/18:2",pgColumn);
      int pgSumColumn = 3;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PG","Overall",pgSumColumn);
      int pe36_2ColumnNa = 4;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PE_Na","18:0/18:2-18:1/18:1",pe36_2ColumnNa);
      int pe36_2ColumnH = 5;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PE_H","18:0/18:2-18:1/18:1",pe36_2ColumnH);
      int pe36_2Column = 6;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PE","18:0/18:2-18:1/18:1",pe36_2Column);
      int pe36_4ColumnNa = 7;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PE_Na","16:0/20:4-18:2/18:2",pe36_4ColumnNa);
      int pe36_4ColumnH = 8;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PE_H","16:0/20:4-18:2/18:2",pe36_4ColumnH);
      int pe36_4Column = 9;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PE","16:0/20:4-18:2/18:2",pe36_4Column);
      int peSumColumn = 10;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PE","Overall",peSumColumn);
      int pc32_0Column_Na = 11;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_Na","14:0/18:0-16:0/16:0",pc32_0Column_Na);
      int pc32_0Column_H = 12;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_H","14:0/18:0-16:0/16:0",pc32_0Column_H);
      int pc32_0Column = 13;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC","14:0/18:0-16:0/16:0",pc32_0Column);
      int pc34_0Column_Na = 14;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_Na","16:0/18:0-17:0/17:0",pc34_0Column_Na);
      int pc34_0Column_H = 15;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_H","16:0/18:0-17:0/17:0",pc34_0Column_H);
      int pc34_0Column = 16;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC","16:0/18:0-17:0/17:0",pc34_0Column);
      int pc36_2Column_Na = 17;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_Na","18:0/18:2-18:1/18:1",pc36_2Column_Na);
      int pc36_2Column_H = 18;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_H","18:0/18:2-18:1/18:1",pc36_2Column_H);
      int pc36_2Column = 19;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC","18:0/18:2-18:1/18:1",pc36_2Column);
      int pcSumColumn = 20;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PC","Overall",pcSumColumn);
      int tgColumn_Na = 21;
      setHeaderCells(sheet,row,secondRow,headerStyle,"TG_Na","16:0/18:0/16:0-19:0/12:0/19:0",tgColumn_Na);
      int tgColumn_NH4 = 22;
      setHeaderCells(sheet,row,secondRow,headerStyle,"TG_NH4","16:0/18:0/16:0-19:0/12:0/19:0",tgColumn_NH4);
      int tgColumn = 23;
      setHeaderCells(sheet,row,secondRow,headerStyle,"TG","16:0/18:0/16:0-19:0/12:0/19:0",tgColumn);
      int tgSumColumn = 24;
      setHeaderCells(sheet,row,secondRow,headerStyle,"TG","Overall",tgSumColumn);
      int psColumn = 25;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PS","18:0/18:2-18:1/18:1",psColumn);
      int psSumColumn = 26;
      setHeaderCells(sheet,row,secondRow,headerStyle,"PS","Overall",psSumColumn);

      
      //negative columns
//      int pgColumn = 2;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PG","18:1/18:1-18:0/18:2",pgColumn);
//      int pgSumColumn = 3;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PG","Overall",pgSumColumn);
//      int pe36_2Column = 4;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PE","18:0/18:2-18:1/18:1",pe36_2Column);
//      int pe36_4Column = 5;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PE","16:0/20:4-18:2/18:2",pe36_4Column);
//      int peSumColumn = 6;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PE","Overall",peSumColumn);
//      int pc32_0Column_CH3 = 7;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_-CH3","14:0/18:0-16:0/16:0",pc32_0Column_CH3);
//      int pc32_0Column_HCOO = 8;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_HCOO","14:0/18:0-16:0/16:0",pc32_0Column_HCOO);
//      int pc32_0Column = 9;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC","14:0/18:0-16:0/16:0",pc32_0Column);
//      int pc34_0Column_CH3 = 10;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_-CH3","16:0/18:0-17:0/17:0",pc34_0Column_CH3);
//      int pc34_0Column_HCOO = 11;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_HCOO","16:0/18:0-17:0/17:0",pc34_0Column_HCOO);
//      int pc34_0Column = 12;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC","16:0/18:0-17:0/17:0",pc34_0Column);
//      int pc36_2Column_CH3 = 13;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_-CH3","18:0/18:2-18:1/18:1",pc36_2Column_CH3);
//      int pc36_2Column_HCOO = 14;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC_HCOO","18:0/18:2-18:1/18:1",pc36_2Column_HCOO);
//      int pc36_2Column = 15;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC","18:0/18:2-18:1/18:1",pc36_2Column);
//      int pcSumColumn = 16;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PC","Overall",pcSumColumn);
//      int psColumn = 17;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PS","18:0/18:2-18:1/18:1",psColumn);
//      int psSumColumn = 18;
//      setHeaderCells(sheet,row,secondRow,headerStyle,"PS","Overall",psSumColumn);

      
      rowCount++;
      
      Vector<String> files = new Vector<String>();
      files.add(baseDir+"002_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"003_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"004_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"005_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"006_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"034_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"035_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"036_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"037_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"038_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-001_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-002_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-003_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-004_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-005_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-001_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-002_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-003_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-004_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-005_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"002_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"003_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"004_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"005_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"006_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"034_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"035_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"036_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"037_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"038_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_05_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_01n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_02n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_03n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_04n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_1_05n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_1a_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_1a_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_1a_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_1a_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_1a_05_Ex3_pos.xlsx");

      
      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,pe36_2ColumnNa,pe36_2ColumnH,pe36_2Column,pe36_4ColumnNa,pe36_4ColumnH,
          pe36_4Column,peSumColumn,pc32_0Column_Na,pc32_0Column_H,pc32_0Column,pc34_0Column_Na,pc34_0Column_H,pc34_0Column,pc36_2Column_Na,
          pc36_2Column_H,pc36_2Column,pcSumColumn,tgColumn_Na,tgColumn_NH4,tgColumn,tgSumColumn,psColumn,psSumColumn,true);
//      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,-1,-1,pe36_2Column,-1,-1,
//          pe36_4Column,peSumColumn,pc32_0Column_CH3,pc32_0Column_HCOO,pc32_0Column,pc34_0Column_CH3,pc34_0Column_HCOO,pc34_0Column,pc36_2Column_CH3,
//          pc36_2Column_HCOO,pc36_2Column,pcSumColumn,-1,-1,-1,-1,psColumn,psSumColumn,false);

      rowCount++;
      files = new Vector<String>();
      files.add(baseDir+"008_Ex3-2_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"009_Ex3-2_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"010_Ex3-2_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"011_Ex3-2_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"012_Ex3-2_Orbitrap_CID_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"040_Ex3-2_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"041_Ex3-2_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"042_Ex3-2_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"043_Ex3-2_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"044_Ex3-2_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"008_Ex3-2_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"009_Ex3-2_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"010_Ex3-2_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"011_Ex3-2_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"012_Ex3-2_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-007_QTrap_Ex3.2_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-008_QTrap_Ex3.2_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-009_QTrap_Ex3.2_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-010_QTrap_Ex3.2_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-011_QTrap_Ex3.2_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-007_QTrap_Ex3.2_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-008_QTrap_Ex3.2_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-009_QTrap_Ex3.2_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-010_QTrap_Ex3.2_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-011_QTrap_Ex3.2_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"040_Ex3-2_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"041_Ex3-2_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"042_Ex3-2_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"043_Ex3-2_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"044_Ex3-2_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_05_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_01n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_02n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_03n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_04n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_2_05n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_2a_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_2a_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_2a_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_2a_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_2a_05_Ex3_pos.xlsx");
      
      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,pe36_2ColumnNa,pe36_2ColumnH,pe36_2Column,pe36_4ColumnNa,pe36_4ColumnH,
          pe36_4Column,peSumColumn,pc32_0Column_Na,pc32_0Column_H,pc32_0Column,pc34_0Column_Na,pc34_0Column_H,pc34_0Column,pc36_2Column_Na,
          pc36_2Column_H,pc36_2Column,pcSumColumn,tgColumn_Na,tgColumn_NH4,tgColumn,tgSumColumn,psColumn,psSumColumn,true);
//      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,-1,-1,pe36_2Column,-1,-1,
//          pe36_4Column,peSumColumn,pc32_0Column_CH3,pc32_0Column_HCOO,pc32_0Column,pc34_0Column_CH3,pc34_0Column_HCOO,pc34_0Column,pc36_2Column_CH3,
//          pc36_2Column_HCOO,pc36_2Column,pcSumColumn,-1,-1,-1,-1,psColumn,psSumColumn,false);

      rowCount++;
      files = new Vector<String>();
      files.add(baseDir+"014_Ex3-3_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"015_Ex3-3_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"016_Ex3-3_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"017_Ex3-3_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"018_Ex3-3_Orbitrap_CID_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"046_Ex3-3_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"047_Ex3-3_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"048_Ex3-3_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"049_Ex3-3_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"050_Ex3-3_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"014_Ex3-3_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"015_Ex3-3_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"016_Ex3-3_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"017_Ex3-3_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"018_Ex3-3_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-013_QTrap_Ex3.3_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-014_QTrap_Ex3.3_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-015_QTrap_Ex3.3_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-016_QTrap_Ex3.3_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-017_QTrap_Ex3.3_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-013_QTrap_Ex3.3_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-014_QTrap_Ex3.3_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-015_QTrap_Ex3.3_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-016_QTrap_Ex3.3_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-017_QTrap_Ex3.3_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"046_Ex3-3_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"047_Ex3-3_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"048_Ex3-3_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"049_Ex3-3_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"050_Ex3-3_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_05_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_01n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_02n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_03n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_04n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_3_05n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_3a_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_3a_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_3a_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_3a_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151125_GNR_LDA_Exp3_3a_05_Ex3_pos.xlsx");

      
      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,pe36_2ColumnNa,pe36_2ColumnH,pe36_2Column,pe36_4ColumnNa,pe36_4ColumnH,
          pe36_4Column,peSumColumn,pc32_0Column_Na,pc32_0Column_H,pc32_0Column,pc34_0Column_Na,pc34_0Column_H,pc34_0Column,pc36_2Column_Na,
          pc36_2Column_H,pc36_2Column,pcSumColumn,tgColumn_Na,tgColumn_NH4,tgColumn,tgSumColumn,psColumn,psSumColumn,true);
//      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,-1,-1,pe36_2Column,-1,-1,
//          pe36_4Column,peSumColumn,pc32_0Column_CH3,pc32_0Column_HCOO,pc32_0Column,pc34_0Column_CH3,pc34_0Column_HCOO,pc34_0Column,pc36_2Column_CH3,
//          pc36_2Column_HCOO,pc36_2Column,pcSumColumn,-1,-1,-1,-1,psColumn,psSumColumn,false);

      rowCount++;
      files = new Vector<String>();
      files.add(baseDir+"020_Ex3-4_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"021_Ex3-4_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"022_Ex3-4_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"023_Ex3-4_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"024_Ex3-4_Orbitrap_CID_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"052_Ex3-4_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"053_Ex3-4_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"054_Ex3-4_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"055_Ex3-4_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"056_Ex3-4_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-019_QTrap_Ex3.4_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-020_QTrap_Ex3.4_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-021_QTrap_Ex3.4_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-022_QTrap_Ex3.4_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-023_QTrap_Ex3.4_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-019_QTrap_Ex3.4_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-020_QTrap_Ex3.4_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-021_QTrap_Ex3.4_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-022_QTrap_Ex3.4_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-023_QTrap_Ex3.4_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"020_Ex3-4_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"021_Ex3-4_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"022_Ex3-4_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"023_Ex3-4_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"024_Ex3-4_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"052_Ex3-4_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"053_Ex3-4_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"054_Ex3-4_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"055_Ex3-4_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"056_Ex3-4_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_05_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_01n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_02n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_03n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_04n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_4_05n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_4_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_4_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_4_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_4_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_4_05_Ex3_pos.xlsx");
      
      
      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,pe36_2ColumnNa,pe36_2ColumnH,pe36_2Column,pe36_4ColumnNa,pe36_4ColumnH,
          pe36_4Column,peSumColumn,pc32_0Column_Na,pc32_0Column_H,pc32_0Column,pc34_0Column_Na,pc34_0Column_H,pc34_0Column,pc36_2Column_Na,
          pc36_2Column_H,pc36_2Column,pcSumColumn,tgColumn_Na,tgColumn_NH4,tgColumn,tgSumColumn,psColumn,psSumColumn,true);
//      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,-1,-1,pe36_2Column,-1,-1,
//          pe36_4Column,peSumColumn,pc32_0Column_CH3,pc32_0Column_HCOO,pc32_0Column,pc34_0Column_CH3,pc34_0Column_HCOO,pc34_0Column,pc36_2Column_CH3,
//          pc36_2Column_HCOO,pc36_2Column,pcSumColumn,-1,-1,-1,-1,psColumn,psSumColumn,false);

      rowCount++;
      files = new Vector<String>();
      files.add(baseDir+"026_Ex3-5_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"027_Ex3-5_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"028_Ex3-5_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"029_Ex3-5_Orbitrap_CID_pos_Ex3_pos.xlsx");
      files.add(baseDir+"030_Ex3-5_Orbitrap_CID_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"058_Ex3-5_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"059_Ex3-5_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"060_Ex3-5_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"061_Ex3-5_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"062_Ex3-5_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-025_QTrap_Ex3.5_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-026_QTrap_Ex3.5_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-027_QTrap_Ex3.5_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-028_QTrap_Ex3.5_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150826_Ex3_QTrap_neg-029_QTrap_Ex3.5_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-025_QTrap_Ex3.5_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-026_QTrap_Ex3.5_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-027_QTrap_Ex3.5_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-028_QTrap_Ex3.5_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"Data20150827_Ex3_QTrap_pos-029_QTrap_Ex3.5_pos_Ex3_pos.xlsx");
//      files.add(baseDir+"026_Ex3-5_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"027_Ex3-5_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"028_Ex3-5_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"029_Ex3-5_Orbitrap_CID_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"030_Ex3-5_Orbitrap_CID_neg_Ex3_neg.xlsx");

//      files.add(baseDir+"058_Ex3-5_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"059_Ex3-5_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"060_Ex3-5_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"061_Ex3-5_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"062_Ex3-5_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_05_Ex3_pos.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_01n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_02n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_03n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_04n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151103_GNR_LDA_Exp3_5_05n_Ex3_neg.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_5_01_Ex3_pos.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_5_02_Ex3_pos.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_5_03_Ex3_pos.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_5_04_Ex3_pos.xlsx");
//      files.add(baseDir+"20151203_LDA_Exp3_5_05_Ex3_pos.xlsx");

      
      
      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,pe36_2ColumnNa,pe36_2ColumnH,pe36_2Column,pe36_4ColumnNa,pe36_4ColumnH,
          pe36_4Column,peSumColumn,pc32_0Column_Na,pc32_0Column_H,pc32_0Column,pc34_0Column_Na,pc34_0Column_H,pc34_0Column,pc36_2Column_Na,
          pc36_2Column_H,pc36_2Column,pcSumColumn,tgColumn_Na,tgColumn_NH4,tgColumn,tgSumColumn,psColumn,psSumColumn,true);
//      rowCount = evaluateExp3SampleBatch(files,sheet,headerStyle,rightStyle,rowCount,pgColumn,pgSumColumn,-1,-1,pe36_2Column,-1,-1,
//          pe36_4Column,peSumColumn,pc32_0Column_CH3,pc32_0Column_HCOO,pc32_0Column,pc34_0Column_CH3,pc34_0Column_HCOO,pc34_0Column,pc36_2Column_CH3,
//          pc36_2Column_HCOO,pc36_2Column,pcSumColumn,-1,-1,-1,-1,psColumn,psSumColumn,false);
      
      rowCount++;
      
      resultWorkbook.write(out);
      resultWorkbook.close();
      out.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void setHeaderCells(Sheet sheet, Row row, Row secondRow, CellStyle headerStyle, String firstValue, String secondValue, int column){
    int width = (19+1)*256;
    sheet.setColumnWidth(column, width);
    Cell cell = row.createCell(column,Cell.CELL_TYPE_STRING);
    cell.setCellStyle(headerStyle);
    cell.setCellValue(firstValue);
    cell = secondRow.createCell(column,Cell.CELL_TYPE_STRING);
    cell.setCellValue(secondValue);

  }
  
  private CellStyle getHeaderStyle(Workbook wb){
    CellStyle arial12style = wb.createCellStyle();
    Font arial12font = wb.createFont();
    arial12font.setBoldweight(Font.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(CellStyle.ALIGN_CENTER);
    return arial12style;
  }

  private CellStyle getLeftHeaderStyle(Workbook wb){
    CellStyle arial12style = getHeaderStyle(wb);
    arial12style.setAlignment(CellStyle.ALIGN_LEFT);
    return arial12style;
  }

  private CellStyle getRightStyle(Workbook wb){
    CellStyle style = wb.createCellStyle();
    style.setAlignment(CellStyle.ALIGN_RIGHT);
    return style;
  }

  private CellStyle getRightIndentStyle(Workbook wb){
    CellStyle style = wb.createCellStyle();
    style.setAlignment(CellStyle.ALIGN_RIGHT);
    short indent = 1;
    style.setIndention(indent);
    return style;
  }
  
  private CellStyle getBoldStyle(Workbook wb){
    CellStyle calibri11style = wb.createCellStyle();
    Font calibri11font = wb.createFont();
    calibri11font.setBoldweight(Font.BOLDWEIGHT_BOLD);
    calibri11font.setFontName("Calibri");
    calibri11font.setFontHeightInPoints((short)11);
    calibri11style.setFont(calibri11font);
    return calibri11style;
  }
  
  private CellStyle getPercentageStyle(Workbook wb){
    CellStyle style = wb.createCellStyle();
    style.setDataFormat(wb.createDataFormat().getFormat("0%"));
    return style;
  }

  
  private CellStyle getCenterStyle(Workbook wb){
    CellStyle style = wb.createCellStyle();
    style.setAlignment(CellStyle.ALIGN_CENTER);
    return style;
  }
  
  // !!!!!!!!!!!!!!!!! ATTENTION: for PE 36:4 the retention time has to be set !!!!!!!!!!!!!!!!!!!
  private int evaluateExp3SampleBatch(Vector<String> files, Sheet sheet, CellStyle headerStyle, CellStyle rightStyle, int count, int pgColumn, int pgSumColumn, int pe36_2ColumnNa, int pe36_2ColumnH,
      int pe36_2Column, int pe36_4ColumnNa, int pe36_4ColumnH, int pe36_4Column, int peSumColumn, int pc32_0ColumnNa, int pc32_0ColumnH,
      int pc32_0Column, int pc34_0ColumnNa, int pc34_0ColumnH, int pc34_0Column, int pc36_2ColumnNa, int pc36_2ColumnH,
      int pc36_2Column, int pcSumColumn, int tgColumnNa, int tgColumnNH4, int tgColumn, int tgSumColumn, int psColumn, int psSumColumn, boolean positive) throws Exception{
    int rowCount=count;
    
    Vector<Float> pgRatios = new Vector<Float>();
    Vector<Float> pe362NaRatios = new Vector<Float>();
    Vector<Float> pe362HRatios = new Vector<Float>();
    Vector<Float> pe362Ratios = new Vector<Float>();
    Vector<Float> pe364NaRatios = new Vector<Float>();
    Vector<Float> pe364HRatios = new Vector<Float>();
    Vector<Float> pe364Ratios = new Vector<Float>();
    Vector<Float> peSumRatios = new Vector<Float>();
    Vector<Float> pc320NaRatios = new Vector<Float>();
    Vector<Float> pc320HRatios = new Vector<Float>();
    Vector<Float> pc320Ratios = new Vector<Float>();
    Vector<Float> pc340NaRatios = new Vector<Float>();
    Vector<Float> pc340HRatios = new Vector<Float>();
    Vector<Float> pc340Ratios = new Vector<Float>();
    Vector<Float> pc362NaRatios = new Vector<Float>();
    Vector<Float> pc362HRatios = new Vector<Float>();
    Vector<Float> pc362Ratios = new Vector<Float>();
    Vector<Float> pcSumRatios = new Vector<Float>();
    Vector<Float> tgNaRatios = new Vector<Float>();
    Vector<Float> tgNH4Ratios = new Vector<Float>();
    Vector<Float> tgRatios = new Vector<Float>();
    Vector<Float> psRatios = new Vector<Float>();
    
    int detectedPGs = 0;
    int detectablePGs = 1;
    int detected362HPEs = 0;
    int detected362NaPEs = 0;
    int detected364HPEs = 0;
    int detected364NaPEs = 0;
    int detected362PEs = 0;
    int detected364PEs = 0;
    int detectablePEs = 2;   
    int detected320NaPCs = 0;
    int detected320HPCs = 0;
    int detected320PCs = 0;
    int detected340NaPCs = 0;
    int detected340HPCs = 0;
    int detected340PCs = 0;
    int detected362NaPCs = 0;
    int detected362HPCs = 0;
    int detected362PCs = 0;
    int detectablePCs = 3;
    int detectedNaTGs = 0;
    int detectedNH4TGs = 0;
    int detectedTGs = 0;
    int detectableTGs = 1;
    int detectedPSs = 0;
    int detectablePSs = 1;
    
    for (String file : files){
      Row row = sheet.createRow(rowCount);
      String fileName = file.substring(file.lastIndexOf("\\")+1);
      Hashtable<String,Vector<LipidParameterSet>> result = LDAResultReader.readResultFile(file, new Hashtable<String,Boolean>()).getIdentifications();
      Cell cell = row.createCell(0);
      cell.setCellType(Cell.CELL_TYPE_STRING);
      cell.setCellValue(fileName);
      
      if (result.containsKey("PG")){
        Vector<LipidParameterSet> sets = result.get("PG");
        Hashtable<String,String> foundPGSpecies= new Hashtable<String,String>();
        float firstArea = 0f;
        float secondArea = 0f;
        boolean firstPGFound = false;
        boolean secondPGFound = false;
        for (LipidParameterSet set: sets){
          if (positive){
            if (!set.getModificationName().equalsIgnoreCase("Na") || !(set instanceof LipidomicsMSnSet)) continue;
          } else {
            if (!set.getModificationName().equalsIgnoreCase("-H") || !(set instanceof LipidomicsMSnSet)) continue;
          }
          float[] areas = separateMSnPeak((LipidomicsMSnSet)set, new String[]{"18:1"}, new String[]{"18:0","18:2"});
          if (areas[0]>0 && areas[1]>0){
//for QTOF         
          ////          if (areas[0]>0 || areas[1]>0){
            firstArea += areas[0];
            secondArea += areas[1];
          }
          if (areas[0]>0) firstPGFound = true;
          if (areas[1]>0) secondPGFound = true;          
        }
        if (firstArea>0 && secondArea>0){
          float ratio = firstArea/secondArea;
          cell = row.createCell(pgColumn);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pgRatios.add(ratio);
        }
        if (firstPGFound) detectedPGs++;
        if (secondPGFound) detectedPGs++;
      }
      if (result.containsKey("PE")){
        Vector<LipidParameterSet> sets = result.get("PE");
        float first362AreaNa = 0f;
        float second362AreaNa = 0f;
        float first362AreaH = 0f;
        float second362AreaH = 0f;
        float first362Area = 0f;
        float second362Area = 0f;
        float first364AreaNa = 0f;
        float second364AreaNa = 0f;
        float first364AreaH = 0f;
        float second364AreaH = 0f;
        float first364Area = 0f;
        float second364Area = 0f;
        
        boolean first362NaFound = false;
        boolean second362NaFound = false;
        boolean first362HFound = false;
        boolean second362HFound = false;
        boolean first362Found = false;
        boolean second362Found = false;
        boolean first364NaFound = false;
        boolean second364NaFound = false;
        boolean first364HFound = false;
        boolean second364HFound = false;
        boolean first364Found = false;
        boolean second364Found = false;
        
        for (LipidParameterSet set: sets){
          if (!(set instanceof LipidomicsMSnSet)) continue;
          if (positive){
            if (!set.getModificationName().equalsIgnoreCase("Na") && !set.getModificationName().equalsIgnoreCase("H")) continue;
          }else{
            if (!set.getModificationName().equalsIgnoreCase("-H")) continue;
          }
          if (set.getNameStringWithoutRt().equalsIgnoreCase("36:2")){
            float[] areas = separateMSnPeak((LipidomicsMSnSet)set, new String[]{"18:0","18:2"}, new String[]{"18:1"});
            if (areas[0]>0 && areas[1]>0){
//for QTOF
////                        if (areas[0]>0 || areas[1]>0){
              if (set.getModificationName().equalsIgnoreCase("Na")){
                first362AreaNa += areas[0];
                second362AreaNa += areas[1];
              } else if ((positive && set.getModificationName().equalsIgnoreCase("H")) || (!positive && set.getModificationName().equalsIgnoreCase("-H"))){
//              } else if (set.getModificationName().equalsIgnoreCase("-H")){
                first362AreaH += areas[0];
                second362AreaH += areas[1];                
              }
              first362Area += areas[0];
              second362Area += areas[1];
            }
            if (areas[0]>0){
              first362Found = true;
              if (set.getModificationName().equalsIgnoreCase("Na")) first362NaFound = true;
              else if (set.getModificationName().equalsIgnoreCase("H")) first362HFound = true;
            }
            if (areas[1]>0){
              second362Found = true;
              if (set.getModificationName().equalsIgnoreCase("Na")) second362NaFound = true;
              else if (set.getModificationName().equalsIgnoreCase("H")) second362HFound = true;
            }
          } else if (set.getNameStringWithoutRt().equalsIgnoreCase("36:4")){
            float rt = Float.parseFloat(set.getRt());
            //the first one is a special case, where both peaks are united
            if (fileName.equalsIgnoreCase("Data20150826_Ex3_QTrap_neg-003_QTrap_Ex3.1_neg_Ex3_neg.xlsx")||
                fileName.equalsIgnoreCase("Data20150826_Ex3_QTrap_neg-015_QTrap_Ex3.3_neg_Ex3_neg.xlsx")||
                fileName.equalsIgnoreCase("Data20150826_Ex3_QTrap_neg-025_QTrap_Ex3.5_neg_Ex3_neg.xlsx")||
                fileName.equalsIgnoreCase("Data20150826_Ex3_QTrap_neg-029_QTrap_Ex3.5_neg_Ex3_neg.xlsx")){
              float[] areas = separateMSnPeak((LipidomicsMSnSet)set, new String[]{"16:0","20:4"}, new String[]{"18:2"});
              if ((positive && set.getModificationName().equalsIgnoreCase("H")) || (!positive && set.getModificationName().equalsIgnoreCase("-H"))){
                first364Area += areas[0];
                second364Area += areas[1];
                if (areas[0]>0) first364Found = true;
                if (areas[1]>0) second364Found = true;
              }
            } else if (22.6f<rt && rt<24.3f){
              if (set.getModificationName().equalsIgnoreCase("Na")){
                second364AreaNa = set.Area;
                second364NaFound = true;
              } else if (set.getModificationName().equalsIgnoreCase("H")){
                second364AreaH = set.Area;
                second364HFound = true;
              }
              second364Area += set.Area;
              second364Found = true;
            } else if (24.3f<rt && rt<25.8f){
              if (set.getModificationName().equalsIgnoreCase("Na")){
                first364AreaNa = set.Area;
                first364NaFound = true;
              }else if (set.getModificationName().equalsIgnoreCase("H")){
                first364AreaH = set.Area;
                first364HFound = true;
              } else if (set.getModificationName().equalsIgnoreCase("-H")){
                first364AreaH = set.Area;
                first364HFound = true;
              }
              first364Area += set.Area;
              first364Found = true;
            } 
          }
        }
        if (positive && first362AreaNa>0 && second362AreaNa>0){
          float ratio = first362AreaNa/second362AreaNa;
          cell = row.createCell(pe36_2ColumnNa);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pe362NaRatios.add(ratio);
        } 
        if (positive && first362NaFound) detected362NaPEs += 1;
        if (positive && second362NaFound) detected362NaPEs += 1;
        if (positive && first362AreaH>0 && second362AreaH>0){
          float ratio = first362AreaH/second362AreaH;
          cell = row.createCell(pe36_2ColumnH);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pe362HRatios.add(ratio);
        }
        if (positive && first362HFound) detected362HPEs += 1;
        if (positive && second362HFound) detected362HPEs += 1;
        if (first362Area>0 && second362Area>0){
          float ratio = first362Area/second362Area;
          cell = row.createCell(pe36_2Column);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pe362Ratios.add(ratio);
          peSumRatios.add(ratio);
        } 
        if (first362Found) detected362PEs += 1;
        if (second362Found) detected362PEs += 1;
        if (positive && first364AreaNa>0 && second364AreaNa>0){
          float ratio = first364AreaNa/second364AreaNa;
          cell = row.createCell(pe36_4ColumnNa);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pe364NaRatios.add(ratio);
        }
        if (positive && first364NaFound) detected364NaPEs += 1;
        if (positive && second364NaFound) detected364NaPEs += 1;
        if (positive && first364AreaH>0 && second364AreaH>0){
          float ratio = first364AreaH/second364AreaH;
          cell = row.createCell(pe36_4ColumnH);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pe364HRatios.add(ratio);
        } 
        if (positive && first364HFound) detected364HPEs += 1;
        if (positive && second364HFound) detected364HPEs += 1;
        if (first364Area>0 && second364Area>0){
          float ratio = first364Area/second364Area;
          cell = row.createCell(pe36_4Column);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pe364Ratios.add(ratio);
          peSumRatios.add(ratio);
        }
        if (first364Found) detected364PEs += 1;
        if (second364Found) detected364PEs += 1;

      }
      if (result.containsKey("PC")){
        Vector<LipidParameterSet> sets = result.get("PC");
        float first320AreaNa = 0f;
        float second320AreaNa = 0f;
        float first320AreaH = 0f;
        float second320AreaH = 0f;
        float first320Area = 0f;
        float second320Area = 0f;
        float first340AreaNa = 0f;
        float second340AreaNa = 0f;
        float first340AreaH = 0f;
        float second340AreaH = 0f;
        float first340Area = 0f;
        float second340Area = 0f;
        float first362AreaNa = 0f;
        float second362AreaNa = 0f;
        float first362AreaH = 0f;
        float second362AreaH = 0f;
        float first362Area = 0f;
        float second362Area = 0f;
        
        boolean first320NaFound = false;
        boolean second320NaFound = false;
        boolean first320HFound = false;
        boolean second320HFound = false;
        boolean first320Found = false;
        boolean second320Found = false;
        boolean first340NaFound = false;
        boolean second340NaFound = false;
        boolean first340HFound = false;
        boolean second340HFound = false;
        boolean first340Found = false;
        boolean second340Found = false;
        boolean first362NaFound = false;
        boolean second362NaFound = false;
        boolean first362HFound = false;
        boolean second362HFound = false;
        boolean first362Found = false;
        boolean second362Found = false;

        for (LipidParameterSet set: sets){
          if (!(set instanceof LipidomicsMSnSet)) continue;
          if (positive){
            if (!set.getModificationName().equalsIgnoreCase("Na") && !set.getModificationName().equalsIgnoreCase("H")) continue;
          } else {
            if (!set.getModificationName().equalsIgnoreCase("HCOO") && !set.getModificationName().equalsIgnoreCase("-CH3")) continue;
          }
          if (set.getNameStringWithoutRt().equalsIgnoreCase("32:0")){
            float[] areas = separateMSnPeak((LipidomicsMSnSet)set, new String[]{"14:0","18:0"}, new String[]{"16:0"});
            if (areas[0]>0 && areas[1]>0){
              if (set.getModificationName().equalsIgnoreCase("Na") || set.getModificationName().equalsIgnoreCase("-CH3")){
//              if (set.getModificationName().equalsIgnoreCase("-CH3")){
                first320AreaNa += areas[0];
                second320AreaNa += areas[1];
              } else if (set.getModificationName().equalsIgnoreCase("H")||set.getModificationName().equalsIgnoreCase("HCOO")){
//              } else if (set.getModificationName().equalsIgnoreCase("HCOO")){
                first320AreaH += areas[0];
                second320AreaH += areas[1];                
              }
              first320Area += areas[0];
              second320Area += areas[1];
            }
            if (areas[0]>0){
              first320Found = true;
              if (set.getModificationName().equalsIgnoreCase("Na") || set.getModificationName().equalsIgnoreCase("-CH3")) first320NaFound = true;
              else if (set.getModificationName().equalsIgnoreCase("H") || set.getModificationName().equalsIgnoreCase("HCOO")) first320HFound = true;
            }
            if (areas[1]>0){
              second320Found = true;
              if (set.getModificationName().equalsIgnoreCase("Na") || set.getModificationName().equalsIgnoreCase("-CH3")) second320NaFound = true;
              else if (set.getModificationName().equalsIgnoreCase("H") || set.getModificationName().equalsIgnoreCase("HCOO")) second320HFound = true;
            }
          } else if (set.getNameStringWithoutRt().equalsIgnoreCase("34:0")){
            float[] areas = separateMSnPeak((LipidomicsMSnSet)set, new String[]{"16:0","18:0"}, new String[]{"17:0"});
            if (areas[0]>0 && areas[1]>0){
            if (set.getModificationName().equalsIgnoreCase("Na")||set.getModificationName().equalsIgnoreCase("-CH3")){
//            if (set.getModificationName().equalsIgnoreCase("-CH3")){
                first340AreaNa += areas[0];
                second340AreaNa += areas[1];
              } else if (set.getModificationName().equalsIgnoreCase("H")||set.getModificationName().equalsIgnoreCase("HCOO")){
//              } else if (set.getModificationName().equalsIgnoreCase("HCOO")){
                first340AreaH += areas[0];
                second340AreaH += areas[1];                
              }
              first340Area += areas[0];
              second340Area += areas[1];
            }
            if (areas[0]>0){
              first340Found = true;
              if (set.getModificationName().equalsIgnoreCase("Na") || set.getModificationName().equalsIgnoreCase("-CH3")) first340NaFound = true;
              else if (set.getModificationName().equalsIgnoreCase("H") || set.getModificationName().equalsIgnoreCase("HCOO")) first340HFound = true;
            }
            if (areas[1]>0){
              second340Found = true;
              if (set.getModificationName().equalsIgnoreCase("Na") || set.getModificationName().equalsIgnoreCase("-CH3")) second340NaFound = true;
              else if (set.getModificationName().equalsIgnoreCase("H") || set.getModificationName().equalsIgnoreCase("HCOO")) second340HFound = true;
            }
          } else if (set.getNameStringWithoutRt().equalsIgnoreCase("36:2")){
            float[] areas = separateMSnPeak((LipidomicsMSnSet)set, new String[]{"18:0","18:2"}, new String[]{"18:1"});
            if (areas[0]>0 && areas[1]>0){
// for QTOF
////            if (areas[0]>0 || areas[1]>0){ 
              if (set.getModificationName().equalsIgnoreCase("Na")||set.getModificationName().equalsIgnoreCase("-CH3")){
//            if (set.getModificationName().equalsIgnoreCase("-CH3")){
                first362AreaNa += areas[0];
                second362AreaNa += areas[1];
              } else if (set.getModificationName().equalsIgnoreCase("H")||set.getModificationName().equalsIgnoreCase("HCOO")){
//              } else if (set.getModificationName().equalsIgnoreCase("HCOO")){
                first362AreaH += areas[0];
                second362AreaH += areas[1];                
              }
              first362Area += areas[0];
              second362Area += areas[1];
            }
            if (areas[0]>0){
              first362Found = true;
              if (set.getModificationName().equalsIgnoreCase("Na") || set.getModificationName().equalsIgnoreCase("-CH3")) first362NaFound = true;
              else if (set.getModificationName().equalsIgnoreCase("H") || set.getModificationName().equalsIgnoreCase("HCOO")) first362HFound = true;
            }
            if (areas[1]>0){
              second362Found = true;
              if (set.getModificationName().equalsIgnoreCase("Na") || set.getModificationName().equalsIgnoreCase("-CH3")) second362NaFound = true;
              else if (set.getModificationName().equalsIgnoreCase("H") || set.getModificationName().equalsIgnoreCase("HCOO")) second362HFound = true;
            }
          }
        }
        if (first320AreaNa>0 && second320AreaNa>0){
          float ratio = first320AreaNa/second320AreaNa;
          cell = row.createCell(pc32_0ColumnNa);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pc320NaRatios.add(ratio);
        }
        if (first320NaFound) detected320NaPCs += 1;
        if (second320NaFound) detected320NaPCs += 1;

        if (first320AreaH>0 && second320AreaH>0){
          float ratio = first320AreaH/second320AreaH;
          cell = row.createCell(pc32_0ColumnH);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pc320HRatios.add(ratio);
        }
        if (first320HFound) detected320HPCs += 1;
        if (second320HFound) detected320HPCs += 1;

        if (first320Area>0 && second320Area>0){
          float ratio = first320Area/second320Area;
          cell = row.createCell(pc32_0Column);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pc320Ratios.add(ratio);
          pcSumRatios.add(ratio);
        }
        if (first320Found) detected320PCs += 1;
        if (second320Found) detected320PCs += 1;

        if (first340AreaNa>0 && second340AreaNa>0){
          float ratio = first340AreaNa/second340AreaNa;
          cell = row.createCell(pc34_0ColumnNa);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pc340NaRatios.add(ratio);
        }
        if (first340NaFound) detected340NaPCs += 1;
        if (second340NaFound) detected340NaPCs += 1;

        if (first340AreaH>0 && second340AreaH>0){
          float ratio = first340AreaH/second340AreaH;
          cell = row.createCell(pc34_0ColumnH);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pc340HRatios.add(ratio);
        }
        if (first340HFound) detected340HPCs += 1;
        if (second340HFound) detected340HPCs += 1;

        if (first340Area>0 && second340Area>0){
          float ratio = first340Area/second340Area;
          cell = row.createCell(pc34_0Column);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pc340Ratios.add(ratio);
          pcSumRatios.add(ratio);
        }
        if (first340Found) detected340PCs += 1;
        if (second340Found) detected340PCs += 1;

        if (first362AreaNa>0 && second362AreaNa>0){
          float ratio = first362AreaNa/second362AreaNa;
          cell = row.createCell(pc36_2ColumnNa);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pc362NaRatios.add(ratio);
        }
        if (first362NaFound) detected362NaPCs += 1;
        if (second362NaFound) detected362NaPCs += 1;

        if (first362AreaH>0 && second362AreaH>0){
          float ratio = first362AreaH/second362AreaH;
          cell = row.createCell(pc36_2ColumnH);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pc362HRatios.add(ratio);
        }
        if (first362HFound) detected362HPCs += 1;
        if (second362HFound) detected362HPCs += 1;

        if (first362Area>0 && second362Area>0){
          float ratio = first362Area/second362Area;
          cell = row.createCell(pc36_2Column);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          pc362Ratios.add(ratio);
          pcSumRatios.add(ratio);
        }
        if (first362Found) detected362PCs += 1;
        if (second362Found) detected362PCs += 1;
      }
      if (result.containsKey("TG")){
        Vector<LipidParameterSet> sets = result.get("TG");
        float firstAreaNa = 0f;
        float secondAreaNa = 0f;
        float firstAreaNH4 = 0f;
        float secondAreaNH4 = 0f;
        float firstArea = 0f;
        float secondArea = 0f;
        
        boolean firstNaFound = false;
        boolean secondNaFound = false;
        boolean firstNH4Found = false;
        boolean secondNH4Found = false;
        boolean firstFound = false;
        boolean secondFound = false;

        for (LipidParameterSet set: sets){
          if (!(set instanceof LipidomicsMSnSet)) continue;
          if (!set.getModificationName().equalsIgnoreCase("Na") && !set.getModificationName().equalsIgnoreCase("NH4")) continue;
          float[] areas = separateMSnPeak((LipidomicsMSnSet)set, new String[]{"16:0","18:0"}, new String[]{"19:0","12:0"});
          if (areas[0]>0 && areas[1]>0){
            if (set.getModificationName().equalsIgnoreCase("Na")){
              firstAreaNa += areas[0];
              secondAreaNa += areas[1];
            } else if (set.getModificationName().equalsIgnoreCase("NH4")){
              firstAreaNH4 += areas[0];
              secondAreaNH4 += areas[1];                
            }
            firstArea += areas[0];
            secondArea += areas[1];
          }
          if (areas[0]>0){
            firstFound = true;
            if (set.getModificationName().equalsIgnoreCase("Na")) firstNaFound = true;
            else if (set.getModificationName().equalsIgnoreCase("NH4")) firstNH4Found = true;
          }
          if (areas[1]>0){
            secondFound = true;
            if (set.getModificationName().equalsIgnoreCase("Na")) secondNaFound = true;
            else if (set.getModificationName().equalsIgnoreCase("NH4")) secondNH4Found = true;
          }
        }
        if (firstAreaNa>0 && secondAreaNa>0){
          float ratio = firstAreaNa/secondAreaNa;
          cell = row.createCell(tgColumnNa);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          tgNaRatios.add(ratio);
        }
        if (firstNaFound) detectedNaTGs++;
        if (secondNaFound) detectedNaTGs++;
        if (firstAreaNH4>0 && secondAreaNH4>0){
          float ratio = firstAreaNH4/secondAreaNH4;
          cell = row.createCell(tgColumnNH4);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          tgNH4Ratios.add(ratio);
        }
        if (firstNH4Found) detectedNH4TGs++;
        if (secondNH4Found) detectedNH4TGs++;
        if (firstArea>0 && secondArea>0){
          float ratio = firstArea/secondArea;
          cell = row.createCell(tgColumn);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          tgRatios.add(ratio);
        }
        if (firstFound) detectedTGs++;
        if (secondFound) detectedTGs++;
      }
      if (result.containsKey("PS")){
        Vector<LipidParameterSet> sets = result.get("PS");
        float firstArea = 0f;
        float secondArea = 0f;
        boolean firstFound = false;
        boolean secondFound = false;
        for (LipidParameterSet set: sets){
          if (!(set instanceof LipidomicsMSnSet)) continue;
          if (positive){
            if (!set.getModificationName().equalsIgnoreCase("H") || !(set instanceof LipidomicsMSnSet)) continue;
          } else {
            if (!set.getModificationName().equalsIgnoreCase("-H") || !(set instanceof LipidomicsMSnSet)) continue;
          }
          float[] areas = separateMSnPeak((LipidomicsMSnSet)set, new String[]{"18:0","18:2"}, new String[]{"18:1"});
          if (areas[0]>0 && areas[1]>0){
            firstArea += areas[0];
            secondArea += areas[1];
          }
          if (areas[0]>0) firstFound = true;
          if (areas[1]>0) secondFound = true;
        }
        if (firstArea>0 && secondArea>0){
          float ratio = firstArea/secondArea;
          cell = row.createCell(psColumn);
          cell.setCellType(Cell.CELL_TYPE_NUMERIC);
          cell.setCellValue(ratio);
          psRatios.add(ratio);
        }
        if (firstFound) detectedPSs++;
        if (secondFound) detectedPSs++;
      }
      rowCount++;
    }
    Row row = sheet.createRow(rowCount);
    Cell cell = row.createCell(0);
    cell.setCellType(Cell.CELL_TYPE_STRING);
    cell.setCellStyle(headerStyle);
    cell.setCellValue("Detected:");
    rowCount++;
    Row meanRow = sheet.createRow(rowCount);
    cell = meanRow.createCell(0);
    cell.setCellType(Cell.CELL_TYPE_STRING);
    cell.setCellStyle(headerStyle);
    cell.setCellValue("Mean:");
    rowCount++;
    Row medRow = sheet.createRow(rowCount);
    cell = medRow.createCell(0);
    cell.setCellType(Cell.CELL_TYPE_STRING);
    cell.setCellStyle(headerStyle);
    cell.setCellValue("Median:");
    rowCount++;
    Row stdevRow = sheet.createRow(rowCount);
    cell = stdevRow.createCell(0);
    cell.setCellType(Cell.CELL_TYPE_STRING);
    cell.setCellStyle(headerStyle);
    cell.setCellValue("StDev:");
    rowCount++;
    Row cvRow = sheet.createRow(rowCount);
    cell = cvRow.createCell(0);
    cell.setCellType(Cell.CELL_TYPE_STRING);
    cell.setCellStyle(headerStyle);
    cell.setCellValue("CV [%]:");
    
    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pgColumn,detectedPGs,2*files.size(),pgRatios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pgSumColumn,detectedPGs,detectablePGs*files.size()*2,pgRatios);  
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pe36_2ColumnNa,detected362NaPEs,2*files.size(),pe362NaRatios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pe36_2ColumnH,detected362HPEs,2*files.size(),pe362HRatios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pe36_2Column,detected362PEs,2*files.size(),pe362Ratios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pe36_4ColumnNa,detected364NaPEs,2*files.size(),pe364NaRatios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pe36_4ColumnH,detected364HPEs,2*files.size(),pe364HRatios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pe36_4Column,detected364PEs,2*files.size(),pe364Ratios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,peSumColumn,detected362PEs+detected364PEs,detectablePEs*files.size()*2,peSumRatios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pc32_0ColumnNa,detected320NaPCs,2*files.size(),pc320NaRatios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pc32_0ColumnH,detected320HPCs,2*files.size(),pc320HRatios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pc32_0Column,detected320PCs,2*files.size(),pc320Ratios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pc34_0ColumnNa,detected340NaPCs,2*files.size(),pc340NaRatios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pc34_0ColumnH,detected340HPCs,2*files.size(),pc340HRatios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pc34_0Column,detected340PCs,2*files.size(),pc340Ratios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pc36_2ColumnNa,detected362NaPCs,2*files.size(),pc362NaRatios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pc36_2ColumnH,detected362HPCs,2*files.size(),pc362HRatios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pc36_2Column,detected362PCs,2*files.size(),pc362Ratios); 
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,pcSumColumn,detected320PCs+detected340PCs+detected362PCs,detectablePCs*files.size()*2,pcSumRatios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,tgColumnNa,detectedNaTGs,2*files.size(),tgNaRatios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,tgColumnNH4,detectedNH4TGs,2*files.size(),tgNH4Ratios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,tgColumn,detectedTGs,2*files.size(),tgRatios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,tgSumColumn,detectedTGs,detectableTGs*files.size()*2,tgRatios);
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,psColumn,detectedPSs,2*files.size(),psRatios);    
    writeMeanAndStDev(rightStyle,row,meanRow,medRow,stdevRow,cvRow,psSumColumn,detectedPSs,detectablePSs*files.size()*2,psRatios);

    rowCount++;
    return rowCount;
  }
  
  private void writeMeanAndStDev(CellStyle rightAlignment, Row countRow, Row meanRow, Row medianRow, Row stDevRow, Row cvRow, int column, int detected, int detectable, Vector<Float> values){
    if (column<0) return;
    if (values.size()==0 && detectable==0) return;
    Cell cell = countRow.createCell(column);
    cell.setCellType(Cell.CELL_TYPE_STRING);
    String percent =  Calculator.FormatNumberToString(((double)detected*100d)/((double)detectable),1);
    cell.setCellValue(String.valueOf(detected)+"/"+String.valueOf(detectable)+" ("+String.valueOf(percent)+"%)");
    cell.setCellStyle(rightAlignment);
    if (values.size()==0) return;
    cell = meanRow.createCell(column);
    cell.setCellType(Cell.CELL_TYPE_NUMERIC);
    float[] inArray = new float[values.size()];
    for (int i=0; i!=values.size(); i++) inArray[i] = values.get(i);
    float mean = Calculator.mean(inArray);
    cell.setCellValue(mean);
    cell = medianRow.createCell(column);
    cell.setCellType(Cell.CELL_TYPE_NUMERIC);
    cell.setCellValue(Calculator.median(values));
    float stdev = Calculator.stddeviation(inArray);
    if (Float.isNaN(stdev) || Float.isInfinite(stdev)) return;
    cell = stDevRow.createCell(column);
    cell.setCellType(Cell.CELL_TYPE_NUMERIC);
    cell.setCellValue(stdev);
    float cv = (stdev*100f)/mean;
    if (Float.isNaN(stdev) || Float.isInfinite(stdev)) return;
    cell = cvRow.createCell(column);
    cell.setCellType(Cell.CELL_TYPE_NUMERIC);
    cell.setCellValue(cv);

  }
  
  private float[] separateMSnPeak(LipidomicsMSnSet msn, String[] fas1, String[] fas2){
    float[] results = new float[2];
    float totalArea = msn.getArea();
    Hashtable<Float,Float> hit1Peaks = new Hashtable<Float,Float>();
    Hashtable<Float,Float> hit2Peaks = new Hashtable<Float,Float>();
    Hashtable<String,Hashtable<String,CgProbe>> frags = msn.getChainFragments();
    for (String fa : fas1){
      if (!frags.containsKey(fa)) continue;
      for (CgProbe peak : frags.get(fa).values()){
        boolean found = false;
        Vector<Float> toRemove = new Vector<Float>();
        for (Float mz : hit1Peaks.keySet()){
          if (mz-LipidomicsConstants.getMs2MzTolerance()<peak.Mz && peak.Mz<mz+LipidomicsConstants.getMs2MzTolerance()){
            found = true;
            toRemove.add(mz);
          }
        }
        for (Float remove: toRemove){
          hit1Peaks.remove(remove);
        }
        toRemove = new Vector<Float>();
        for (Float mz : hit2Peaks.keySet()){
          if (mz-LipidomicsConstants.getMs2MzTolerance()<peak.Mz && peak.Mz<mz+LipidomicsConstants.getMs2MzTolerance()){
            found = true;
            toRemove.add(mz);
          }
        }
        for (Float remove: toRemove){
          hit2Peaks.remove(remove);
        }
        if (!found) hit1Peaks.put(peak.Mz, peak.Area);
      }
    }
    for (String fa : fas2){
      if (!frags.containsKey(fa)) continue;
      for (CgProbe peak : frags.get(fa).values()){
        boolean found = false;
        Vector<Float> toRemove = new Vector<Float>();
        for (Float mz : hit1Peaks.keySet()){
          if (mz-LipidomicsConstants.getMs2MzTolerance()<peak.Mz && peak.Mz<mz+LipidomicsConstants.getMs2MzTolerance()){
            found = true;
            toRemove.add(mz);
          }
        }
        for (Float remove: toRemove){
          hit1Peaks.remove(remove);
        }
        toRemove = new Vector<Float>();
        for (Float mz : hit2Peaks.keySet()){
          if (mz-LipidomicsConstants.getMs2MzTolerance()<peak.Mz && peak.Mz<mz+LipidomicsConstants.getMs2MzTolerance()){
            found = true;
            toRemove.add(mz);
          }
        }
        for (Float remove: toRemove){
          hit2Peaks.remove(remove);
        }
        if (!found) hit2Peaks.put(peak.Mz, peak.Area);
      }
    }

    float firstProportion = 0f;
    for (float value : hit1Peaks.values()) firstProportion += value;
    float secondProportion = 0f;
    for (float value : hit2Peaks.values()) secondProportion += value;
    results[0] = totalArea*(firstProportion/(firstProportion+secondProportion));
    results[1] = totalArea*(secondProportion/(firstProportion+secondProportion));
    return results;
  }
  
  private void writeLDAResultsToCompareExcel(){
    try{
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream("D:\\Experiment1\\LipidBLAST\\HepatocytesPositiveEvaluation.xlsxn"));
      Workbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      Hashtable<String,Vector<LipidParameterSet>> results = LDAResultReader.readResultFile("D:\\Experiment1\\LipidBLAST\\014_Ex1_Orbitrap_CID_pos_50_DG_LPC_LPE_PC_PE_PG_PS_TG.xlsx", new Hashtable<String,Boolean>()).getIdentifications();
      for (String lClass : results.keySet()){
        Sheet sheet = resultWorkbook.createSheet(lClass);
        int rowCount = 0;
        Row row =  sheet.createRow(rowCount);
        int ms1Column = 0;
        int rtColumn = 1;
        int ldaIdentification = 2;
        rowCount++;
        Cell cell = row.createCell(ms1Column);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Name");
        cell = row.createCell(rtColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("RTime");
        cell = row.createCell(ldaIdentification);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA");
        Hashtable<String,Hashtable<String,String>> usedMs1 = new Hashtable<String,Hashtable<String,String>>(); 
        Vector<LipidParameterSet> identifications = results.get(lClass);
        for (LipidParameterSet set : identifications){
          if (!(set instanceof LipidomicsMSnSet)) continue;
          LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
          for (Object msnNames : msn.getMSnIdentificationNames()){    
            row =  sheet.createRow(rowCount);
            rowCount++;
            String nameString = "";
            if (msnNames instanceof Vector){
              for (String name : (Vector<String>)msnNames){
                nameString += name+";";
              }
              nameString = nameString.substring(0,nameString.length()-1);
            }else{
              nameString = (String)msnNames;
            }
            Hashtable<String,String> rts = new Hashtable<String,String>();
            if (usedMs1.containsKey(set.getNameStringWithoutRt())) rts = usedMs1.get(set.getNameStringWithoutRt());
            else{
              cell = row.createCell(ms1Column);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(set.getNameStringWithoutRt());
            }
            if (!rts.containsKey(set.getRt())){
              cell = row.createCell(rtColumn);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(set.getRt());              
            }
            cell = row.createCell(ldaIdentification);
            cell.setCellType(Cell.CELL_TYPE_STRING);
            cell.setCellValue(nameString);
            rts.put(set.getRt(), set.getRt());
            usedMs1.put(set.getNameStringWithoutRt(), rts);
          }
        }
      }
      resultWorkbook.write(out);
      out.close();
    }catch(Exception ex){
      
    }
  }
  
  

  private void compareLDABLASTNaturalProbesNegative(){
    //the first key is the lipid class, the second key the ms1 species name, the third key the structural identification
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>>();
    Hashtable<String,LipidClassInfoVO> lipidClassInfo = new Hashtable<String,LipidClassInfoVO>();
    LinkedHashMap<String,Boolean> adducts = new LinkedHashMap<String,Boolean>();
    
    ////this.getValidOrbitrapCIDSpeciesNegative(lipidClasses,lipidClassInfo,adducts);
    this.getValid4000QTRAPSpeciesNegative(lipidClasses,lipidClassInfo,adducts);

//    String chromFile = "D:\\BiologicalExperiment\\Orbitrap_CID\\negative\\002_liver2-1_Orbitrap_CID_neg.chrom";
//    String quantFile = "D:\\BiologicalExperiment\\massLists\\negative\\negative.xlsx";
//    String ldaFile = "D:\\BiologicalExperiment\\Orbitrap_CID\\negative\\002_liver2-1_Orbitrap_CID_neg_negative.xlsx";
//    String lbFile = "D:\\BiologicalExperiment\\LipidBlast\\negative\\output\\002_liver2-1_Orbitrap_CID_neg_MF10.mgf.tsv";
//    String outputFile = "D:\\BiologicalExperiment\\LipidBlast\\negative\\002_liver2-1_Orbitrap_CID_neg_LB10_comp_generated.xlsx";
    String chromFile = "D:\\BiologicalExperiment\\QTRAP\\negative\\Data20151002_QTrap_Liver-025_QTrap_Liver1-1_neg.chrom";
    String quantFile = "D:\\BiologicalExperiment\\massLists\\negative\\QTRAP\\negative.xlsx";
    String ldaFile = "D:\\BiologicalExperiment\\QTRAP\\negative\\Data20151002_QTrap_Liver-025_QTrap_Liver1-1_neg_negative.xlsx";
    String lbFile = "D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\negative\\output\\Data20151002_QTrap_Liver-025_QTrap_Liver1-1_neg_MF450.mgf.tsv";
    String outputFile = "D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\negative\\Data20151002_QTrap_Liver-025_QTrap_Liver1-1_neg_MF450_comp_generated.xlsx";

    performComparisonOfNaturalProbes(lipidClasses,lipidClassInfo,chromFile,quantFile,ldaFile,lbFile,outputFile);
  }
  
  private void compareLDABLASTNaturalProbesPositive(){
    //the first key is the lipid class, the second key the ms1 species name, the third key the structural identification
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>>();
    Hashtable<String,LipidClassInfoVO> lipidClassInfo = new Hashtable<String,LipidClassInfoVO>();
    LinkedHashMap<String,Boolean> adducts = new LinkedHashMap<String,Boolean>();
    ////this.getValidOrbitrapCIDSpeciesPositive(lipidClasses,lipidClassInfo,adducts);
    this.getValid4000QTRAPSpeciesPositive(lipidClasses,lipidClassInfo,adducts);
    
    
    String chromFile = "D:\\BiologicalExperiment\\QTRAP\\positive\\Data20151002_QTrap_Liver-003_QTrap_Liver1-1_pos.chrom";
    String quantFile = "D:\\BiologicalExperiment\\massLists\\positive\\QTRAP\\positive.xlsx";
    String ldaFile = "D:\\BiologicalExperiment\\QTRAP\\positive\\Data20151002_QTrap_Liver-003_QTrap_Liver1-1_pos_positive.xlsx";
    String lbFile = "D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\positive\\output\\Data20151002_QTrap_Liver-003_QTrap_Liver1-1_pos_MF10.mgf.tsv";
    String outputFile = "D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\positive\\Data20151002_QTrap_Liver-003_QTrap_Liver1-1_pos_LB10_comp_generated.xlsx";
    performComparisonOfNaturalProbes(lipidClasses,lipidClassInfo,chromFile,quantFile,ldaFile,lbFile,outputFile);
  }
  
  private void correctRetentionTimes(double correction, LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses){
    for (LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> eachClass : lipidClasses.values()){
      correctRetentionTimes(correction,eachClass);
    }
  }
  
  private void correctRetentionTimes(double correction, Map<String,LinkedHashMap<String,ReferenceInfoVO>> toCorrect){
    for (LinkedHashMap<String,ReferenceInfoVO> analyte : toCorrect.values()){
      correctRetentionTimes(correction, analyte.values());
    }
  }

  private void correctRetentionTimes(double correction, Collection<ReferenceInfoVO> toCorrect){
    for (ReferenceInfoVO ref : toCorrect){
      for (int i=0; i!=ref.getCorrectRts().length; i++)
        ref.getCorrectRts()[i] = ref.getCorrectRts()[i]+correction;
    }
  }

  
  private void performComparisonOfNaturalProbes(LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses, 
      Hashtable<String,LipidClassInfoVO> lipidClassInfo, String chromFile, String quantFile, String ldaFile, String lbFile, String outputFile){
    try{
      //Vector quantValues = QuantificationThread.parseQuantExcelFile(quantFile,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f);
      Vector quantValues = null;
      String[] chromPaths = StringUtils.getChromFilePaths(chromFile);   
      LipidomicsAnalyzer analyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);

      Hashtable<String,Vector<String>> analyteSequence = (Hashtable<String,Vector<String>>)quantValues.get(1);
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantVOs = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>) quantValues.get(3);
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outputFile));
      Workbook resultWorkbook = new XSSFWorkbook();
      CellStyle fullyCorrectStyle = getFullyCorrectStyle(resultWorkbook);
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      CellStyle leftHeaderStyle = getLeftHeaderStyle(resultWorkbook);
      CellStyle boldStyle = getBoldStyle(resultWorkbook);
      CellStyle notFoundStyle = getNotFoundStyle(resultWorkbook);
      CellStyle ambigStyle = getFACorrectStyle(resultWorkbook);
      CellStyle onlyMS1Style = getMS1FoundStyle(resultWorkbook);
      CellStyle percentageStyle = getPercentageStyle(resultWorkbook);
      leftHeaderStyle.cloneStyleFrom(headerStyle);
      leftHeaderStyle.setAlignment(CellStyle.ALIGN_LEFT);
      Hashtable<String,Vector<LipidParameterSet>> resultsLDA = LDAResultReader.readResultFile(ldaFile, new Hashtable<String,Boolean>()).getIdentifications();
      LipidBLASTParser lBlastParser = new LipidBLASTParser(lbFile);
      lBlastParser.parse();
      Hashtable<String,Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>> resultsLB = lBlastParser.getResults_();
      for (String lClass : lipidClasses.keySet()){
        LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = lipidClasses.get(lClass);
        LipidClassInfoVO classInfo = lipidClassInfo.get(lClass);
        Sheet sheet = resultWorkbook.createSheet(lClass);
        int rowCount = 0;
        Row row =  sheet.createRow(rowCount);
        int classColumn = 0;
        int width = (10+1)*256;
        int ms1Column = 1;
        sheet.setColumnWidth(ms1Column, width);
        int adductColumn = 2;
        width = (12+1)*256;
        int ldaIdentification = 3;
        sheet.setColumnWidth(ldaIdentification, width);
        int lbIdentification = 4;
        width = (12+1)*256;
        sheet.setColumnWidth(lbIdentification, width);
        int ldaMS1CodeColumn = 5;
        width = (12+1)*256;
        sheet.setColumnWidth(ldaMS1CodeColumn, width);
        int lbMS1CodeColumn = 6;
        width = (12+1)*256;
        sheet.setColumnWidth(lbMS1CodeColumn, width);
        int structureColumn = 7;
        int ldaIdentificationMS2 = 8;
        int lbIdentificationMS2 = 9;
        int ldaMS2CodeColumn = 10;
        int lbMS2CodeColumn = 11;
        int amountMS2Columns = 5;
        int rtColumn = 7;
        int rtLDAColumn = 8;
        int rtLBColumn = 9;
        int probLBColumn = 10;
        int commentColumn = 11;
        int ldaMS1IgnoreColumn = 12;
        int lbMS1IgnoreColumn = 13;
        int ldaMS2IgnoreColumn = 14;
        int lbMS2IgnoreColumn = 15;
        
        rowCount++;
        Cell cell = row.createCell(classColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Class");
        cell = row.createCell(ms1Column);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Name");
        cell = row.createCell(adductColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Adduct");
        cell = row.createCell(ldaIdentification);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA");
        cell = row.createCell(lbIdentification);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LipidBlast");
        cell = row.createCell(ldaMS1CodeColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA-Code");
        cell = row.createCell(lbMS1CodeColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LB-Code");
        if (classInfo.getSns()>1){
          width = (22+1)*256;
          sheet.setColumnWidth(structureColumn, width);
          cell = row.createCell(structureColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("Structure");
          width = (12+1)*256;
          sheet.setColumnWidth(ldaIdentificationMS2, width);
          cell = row.createCell(ldaIdentificationMS2);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("LDA");
          width = (22+1)*256;
          sheet.setColumnWidth(lbIdentificationMS2, width);
          cell = row.createCell(lbIdentificationMS2);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("LipidBlast");
          width = (12+1)*256;
          sheet.setColumnWidth(ldaMS2CodeColumn, width);
          cell = row.createCell(ldaMS2CodeColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("LDA-Code");
          width = (12+1)*256;
          sheet.setColumnWidth(lbMS2CodeColumn, width);
          cell = row.createCell(lbMS2CodeColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("LB-Code");
         
          rtColumn += amountMS2Columns;
          rtLDAColumn += amountMS2Columns;
          rtLBColumn += amountMS2Columns;
          probLBColumn += amountMS2Columns;
          commentColumn += amountMS2Columns;
          ldaMS1IgnoreColumn += amountMS2Columns;
          lbMS1IgnoreColumn += amountMS2Columns;
          ldaMS2IgnoreColumn += amountMS2Columns;
          lbMS2IgnoreColumn += amountMS2Columns;
          
          width = (16+1)*256;
          sheet.setColumnWidth(ldaMS2IgnoreColumn, width);
          cell = row.createCell(ldaMS2IgnoreColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("LDAMS2-Ignore");
          width = (16+1)*256;
          sheet.setColumnWidth(lbMS2IgnoreColumn, width);
          cell = row.createCell(lbMS2IgnoreColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("LBMS2-Ignore");
        }
        width = (9+1)*256;
        sheet.setColumnWidth(rtColumn, width);
        cell = row.createCell(rtColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("RT-ideal");
        width = (9+1)*256;
        sheet.setColumnWidth(rtLDAColumn, width);
        cell = row.createCell(rtLDAColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("RT-LDA");
        width = (9+1)*256;
        sheet.setColumnWidth(rtLBColumn, width);
        cell = row.createCell(rtLBColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("RT-LB");
        width = (9+1)*256;
        sheet.setColumnWidth(probLBColumn, width);
        cell = row.createCell(probLBColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LB-Prob");
        width = (26+1)*256;
        sheet.setColumnWidth(commentColumn, width);
        cell = row.createCell(commentColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Comment");      
        width = (16+1)*256;
        sheet.setColumnWidth(ldaMS1IgnoreColumn, width);
        cell = row.createCell(ldaMS1IgnoreColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDAMS1-Ignore");
        width = (16+1)*256;
        sheet.setColumnWidth(lbMS1IgnoreColumn, width);
        cell = row.createCell(lbMS1IgnoreColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LBMS1-Ignore");


        int finalRowNumberAnalytes = 0;
        Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = quantVOs.get(lClass);
        Vector<LipidParameterSet> ldaAnalytes = resultsLDA.get(lClass);
        Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>> lbAnalytes = new Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>();
        if (resultsLB.containsKey(lClass)) lbAnalytes = resultsLB.get(lClass);
        for (String analyte : analyteSequence.get(lClass)){
          if (analyte.startsWith("IS")) continue;
          Hashtable<String,QuantVO> quantsOfAnalyte = quantsOfClass.get(analyte);
          Hashtable<String,LipidBLASTIdentificationVO> lbsOfMod = new Hashtable<String,LipidBLASTIdentificationVO>();
          if (lbAnalytes.containsKey(analyte)) lbsOfMod = lbAnalytes.get(analyte);
          LinkedHashMap<String,ReferenceInfoVO> molecularSpecies = null;
          if (species.containsKey(analyte)) molecularSpecies = species.get(analyte);
          Vector<LdaLBLASTCompareVO> comparesOfOneSpecies = new Vector<LdaLBLASTCompareVO>();
          
//          for (String adduct : quantsOfAnalyte.keySet()){
//          if (!classInfo.getAdducts().containsKey(adduct)) continue;
          for (String adduct : classInfo.getAdducts().keySet()){
            if (!quantsOfAnalyte.containsKey(adduct)) continue;
            Vector<LipidParameterSet> ldaAnalyte = getLDAAnalytes(ldaAnalytes, analyte, adduct, true, true, molecularSpecies.get(adduct).getCorrectRts());
            LipidBLASTIdentificationVO ident = lbsOfMod.get(adduct);
            if (ident!=null || ldaAnalyte.size()>0){
              LdaLBLASTCompareVO ms1Compare = getLdaLBlastSpeciesComparision(analyte,classInfo,molecularSpecies,ldaAnalyte,ident,analyzer,quantsOfAnalyte.get(adduct));
              //TODO: pay attention for TG this check is not adapted yet
              if (ms1Compare==null) continue;
              if (lClass.equalsIgnoreCase("TG")) checkForUniqueFACombinations(ms1Compare.getMs2Evidence());
              ms1Compare.setAdduct(adduct);
              comparesOfOneSpecies.add(ms1Compare);
            }
          }   
          setIgnoresForSpeciesLevelEvaluation(comparesOfOneSpecies,classInfo.getSns()>1);
          for (LdaLBLASTCompareVO ms1Compare : comparesOfOneSpecies){
            row =  sheet.createRow(rowCount);
            rowCount++;
            cell = row.createCell(classColumn);
            cell.setCellType(Cell.CELL_TYPE_STRING);
            cell.setCellValue(lClass);
            cell = row.createCell(ms1Column);
            cell.setCellType(Cell.CELL_TYPE_STRING);
            if (ms1Compare.isFp()) cell.setCellStyle(notFoundStyle);
            else cell.setCellStyle(fullyCorrectStyle );
            cell.setCellValue(ms1Compare.getCorrectName());
            cell = row.createCell(adductColumn);
            cell.setCellType(Cell.CELL_TYPE_STRING);
            cell.setCellValue(ms1Compare.getAdduct());
              
            CellStyle ldaMS1Style = getStyleBasedOnEvidence(ms1Compare.getLdaIdentCode(),ms1Compare.isFp(),fullyCorrectStyle,
                ambigStyle,onlyMS1Style,notFoundStyle);
            CellStyle lbMS1Style = getStyleBasedOnEvidence(ms1Compare.getLbIdentCode(),ms1Compare.isFp(),fullyCorrectStyle,
                ambigStyle,onlyMS1Style,notFoundStyle);
              
              
              cell = row.createCell(ldaIdentification);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              if (ldaMS1Style!=null) cell.setCellStyle(ldaMS1Style);
              cell.setCellValue(ms1Compare.getLdaName());
              cell = row.createCell(lbIdentification);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              if (lbMS1Style!=null) cell.setCellStyle(lbMS1Style);
              cell.setCellValue(ms1Compare.getLbName());
              cell = row.createCell(ldaMS1CodeColumn);
              cell.setCellType(Cell.CELL_TYPE_NUMERIC);
              if (ldaMS1Style!=null) cell.setCellStyle(ldaMS1Style);
              cell.setCellValue(ms1Compare.getLdaIdentCode());
              cell = row.createCell(lbMS1CodeColumn);
              cell.setCellType(Cell.CELL_TYPE_NUMERIC);
              if (lbMS1Style!=null) cell.setCellStyle(lbMS1Style);
              cell.setCellValue(ms1Compare.getLbIdentCode());
              if (ms1Compare.isIgnoreLDA()){
                cell = row.createCell(ldaMS1IgnoreColumn);
                cell.setCellValue(true);
              }
              if (ms1Compare.isIgnoreLB()){
                cell = row.createCell(lbMS1IgnoreColumn);
                cell.setCellValue(true);
              }
              
              
              if (classInfo.getSns()<2 || ms1Compare.getMs2Evidence().size()==0){
                cell = row.createCell(rtColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(ms1Compare.getRt());
                cell = row.createCell(rtLDAColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(ms1Compare.getLdaRts());
                cell = row.createCell(rtLBColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(ms1Compare.getLbRts());
                if (ms1Compare.getLbIdentCode()>0 || ms1Compare.getLbIdentCode()<-1){
                  cell = row.createCell(probLBColumn);
                  cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(ms1Compare.getLbProb());
                }
                if (ms1Compare.isLdaMs1Only() || ms1Compare.areFasInOtherCombination()){
                  cell = row.createCell(commentColumn);
                  cell.setCellType(Cell.CELL_TYPE_STRING);
                  if (ms1Compare.isLdaMs1Only())
                    cell.setCellValue("not counted: only MS1 identification");
                  else
                    cell.setCellValue("not counted: only possible as additional combination");                  
                }

              }
              
              for (int i=0; i!=ms1Compare.getMs2Evidence().size(); i++){
                LdaLBLASTCompareVO comp = ms1Compare.getMs2Evidence().get(i);
                if (i!=0){
                  row =  sheet.createRow(rowCount);
                  rowCount++;
                }
                cell = row.createCell(structureColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                if (comp.isFp()) cell.setCellStyle(notFoundStyle);
                else cell.setCellStyle(fullyCorrectStyle );
                cell.setCellValue(comp.getCorrectName());

                CellStyle ldaMS2Style = getStyleBasedOnEvidence(comp.getLdaIdentCode(),comp.isFp(),fullyCorrectStyle,
                    ambigStyle,onlyMS1Style,notFoundStyle);
                CellStyle lbMS2Style = getStyleBasedOnEvidence(comp.getLbIdentCode(),comp.isFp(),fullyCorrectStyle,
                    ambigStyle,onlyMS1Style,notFoundStyle);
                cell = row.createCell(ldaIdentificationMS2);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                if (ldaMS1Style!=null) cell.setCellStyle(ldaMS2Style);
                cell.setCellValue(comp.getLdaName());
                cell = row.createCell(lbIdentificationMS2);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                if (lbMS1Style!=null) cell.setCellStyle(lbMS2Style);
                cell.setCellValue(comp.getLbName());
                cell = row.createCell(ldaMS2CodeColumn);
                cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                if (ldaMS2Style!=null) cell.setCellStyle(ldaMS2Style);
                cell.setCellValue(comp.getLdaIdentCode());
                cell = row.createCell(lbMS2CodeColumn);
                cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                if (lbMS1Style!=null) cell.setCellStyle(lbMS2Style);
                cell.setCellValue(comp.getLbIdentCode());
                cell = row.createCell(rtColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(comp.getRt());
                cell = row.createCell(rtLDAColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(comp.getLdaRts());
                cell = row.createCell(rtLBColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(comp.getLbRts());
                if (comp.getLbIdentCode()>0 || comp.getLbIdentCode()<-1){
                  cell = row.createCell(probLBColumn);
                  cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(comp.getLbProb());
                }
                if (comp.isIgnoreLDA()){
                  cell = row.createCell(ldaMS2IgnoreColumn);
                  cell.setCellValue(true);
                }
                if (comp.isIgnoreLB()){
                  cell = row.createCell(lbMS2IgnoreColumn);
                  cell.setCellValue(true);
                }
              }
          }
        }
        finalRowNumberAnalytes = rowCount;

        rowCount++;
        rowCount++;        
        // now starts the summary section of the details tab
        int legendColumn = 0;
        int adductResultColumn = 5;
        int legend2Column = 7;
        int classResultColumn = 10;
        
        String ldaIgnore = "M";
        String lbIgnore = "N";
        if (classInfo.getSns()>1){
          ldaIgnore = "R";
          lbIgnore = "S";          
        }
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Species evaluation");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Species evaluation - adduct insensitive");
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        int totalSumRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified species:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(B2:B"+finalRowNumberAnalytes+",\"<>*_*\",B2:B"+finalRowNumberAnalytes+",\"<>\")");       
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified species:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(B2:B"+finalRowNumberAnalytes+",\"<>*_*\",B2:B"+finalRowNumberAnalytes+",\"<>\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        int ldaCorrectRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified LDA species:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\">0\")");       
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified LDA species:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\">0\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        int lbCorrectRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified LipidBlast species:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\">0\")");       
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified LipidBlast species:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\">0\","+lbIgnore+"2:"+lbIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        int ldaFPRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False positives LDA:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"<>-1\",F2:F"+finalRowNumberAnalytes+",\"<>0\",F2:F"+finalRowNumberAnalytes+",\"<2\")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False positives LDA:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"<>-1\",F2:F"+finalRowNumberAnalytes+",\"<>0\",F2:F"+finalRowNumberAnalytes+",\"<2\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        int lbFPRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False positives LipidBlast:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"<>-1\",G2:G"+finalRowNumberAnalytes+",\"<>0\",G2:G"+finalRowNumberAnalytes+",\"<2\")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False positives LipidBlast:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"<>-1\",G2:G"+finalRowNumberAnalytes+",\"<>0\",G2:G"+finalRowNumberAnalytes+",\"<2\","+lbIgnore+"2:"+lbIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False negatives LDA:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"=-1\")+COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"=-3\")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False negatives LDA:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"=-1\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")+COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"=-3\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False negatives LipidBlast:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"=-1\")+COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"=-3\")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False negatives LipidBlast:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"=-1\","+lbIgnore+"2:"+lbIgnore+finalRowNumberAnalytes+",\"<>TRUE\")+COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"=-3\","+lbIgnore+"2:"+lbIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Sensitivity LDA:");
        cell = row.createCell(adductResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("F"+ldaCorrectRow+"/F"+totalSumRow);
        cell = row.createCell(legend2Column);
        cell.setCellStyle(percentageStyle);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Sensitivity LDA:");
        cell = row.createCell(classResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("K"+ldaCorrectRow+"/K"+totalSumRow);
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Sensitivity LipidBlast");
        cell = row.createCell(adductResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("F"+lbCorrectRow+"/F"+totalSumRow);
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Sensitivity LipidBlast:");
        cell = row.createCell(classResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("K"+lbCorrectRow+"/K"+totalSumRow);
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Positive Predictive Value LDA:");
        cell = row.createCell(adductResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("F"+ldaCorrectRow+"/(F"+ldaCorrectRow+"+F"+ldaFPRow+")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(percentageStyle);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Positive Predictive Value LDA:");
        cell = row.createCell(classResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("K"+ldaCorrectRow+"/(K"+ldaCorrectRow+"+K"+ldaFPRow+")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Positive Predictive Value LipidBlast:");
        cell = row.createCell(adductResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("F"+lbCorrectRow+"/(F"+lbCorrectRow+"+F"+lbFPRow+")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(percentageStyle);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Positive Predictive Value LipidBlast:");
        cell = row.createCell(classResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("K"+lbCorrectRow+"/(K"+lbCorrectRow+"+K"+lbFPRow+")");

        if (classInfo.getSns()>1){
          rowCount++;
          rowCount++;
          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(leftHeaderStyle);
          cell.setCellValue("Molecular species/structure evaluation");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(leftHeaderStyle);
          cell.setCellValue("Molecular species/structure evaluation - adduct insensitive:");
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          totalSumRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified species:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(H2:H"+finalRowNumberAnalytes+",\"<>*_FP\",H2:H"+finalRowNumberAnalytes+",\"<>\",H2:H"+finalRowNumberAnalytes+",\"<>no structure\")");       
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified species:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(H2:H"+finalRowNumberAnalytes+",\"<>*_FP\",H2:H"+finalRowNumberAnalytes+",\"<>\",H2:H"+finalRowNumberAnalytes+",\"<>no structure\",T2:T"+finalRowNumberAnalytes+",\"<>TRUE\")");       
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          ldaCorrectRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified LDA species:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\">0\")");     
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified LDA species:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\">0\",T2:T"+finalRowNumberAnalytes+",\"<>TRUE\")");
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          lbCorrectRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified LipidBlast species:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\">0\")");       
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified LipidBlast species:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\">0\",U2:U"+finalRowNumberAnalytes+",\"<>TRUE\")");
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          ldaFPRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False positives LDA:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"<>-1\",K2:K"+finalRowNumberAnalytes+",\"<>0\",K2:K"+finalRowNumberAnalytes+",\"<2\")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False positives LDA:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"<>-1\",K2:K"+finalRowNumberAnalytes+",\"<>0\",K2:K"+finalRowNumberAnalytes+",\"<2\",T2:T"+finalRowNumberAnalytes+",\"<>TRUE\")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          lbFPRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False positives LipidBlast:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"<>-1\",L2:L"+finalRowNumberAnalytes+",\"<>0\",L2:L"+finalRowNumberAnalytes+",\"<2\")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False positives LipidBlast:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"<>-1\",L2:L"+finalRowNumberAnalytes+",\"<>0\",L2:L"+finalRowNumberAnalytes+",\"<2\",U2:U"+finalRowNumberAnalytes+",\"<>TRUE\")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False negatives LDA:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"=-1\")+COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"=-3\")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False negatives LDA:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"=-1\",T2:T"+finalRowNumberAnalytes+",\"<>TRUE\")+COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"=-3\",T2:T"+finalRowNumberAnalytes+",\"<>TRUE\")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False negatives LipidBlast:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"=-1\")+COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"=-3\")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False negatives LipidBlast:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"=-1\",U2:U"+finalRowNumberAnalytes+",\"<>TRUE\")+COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"=-3\",U2:U"+finalRowNumberAnalytes+",\"<>TRUE\")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Sensitivity LDA:");
          cell = row.createCell(adductResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("F"+ldaCorrectRow+"/F"+totalSumRow);
          cell = row.createCell(legend2Column);
          cell.setCellStyle(percentageStyle);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Sensitivity LDA:");
          cell = row.createCell(classResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("K"+ldaCorrectRow+"/K"+totalSumRow);
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Sensitivity LipidBlast");
          cell = row.createCell(adductResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("F"+lbCorrectRow+"/F"+totalSumRow);
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Sensitivity LipidBlast:");
          cell = row.createCell(classResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("K"+lbCorrectRow+"/K"+totalSumRow);
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Positive Predictive Value LDA:");
          cell = row.createCell(adductResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("F"+ldaCorrectRow+"/(F"+ldaCorrectRow+"+F"+ldaFPRow+")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(percentageStyle);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Positive Predictive Value LDA:");
          cell = row.createCell(classResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("K"+ldaCorrectRow+"/(K"+ldaCorrectRow+"+K"+ldaFPRow+")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Positive Predictive Value LipidBlast:");
          cell = row.createCell(adductResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("F"+lbCorrectRow+"/(F"+lbCorrectRow+"+F"+lbFPRow+")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(percentageStyle);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Positive Predictive Value LipidBlast:");
          cell = row.createCell(classResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("K"+lbCorrectRow+"/(K"+lbCorrectRow+"+K"+lbFPRow+")");
          
        }
        rowCount++;
        rowCount++;
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Legend:");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(fullyCorrectStyle);
        cell.setCellValue("the name is green if the species is correct");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(notFoundStyle);
        cell.setCellValue("the name is red if it is a false positive");
        rowCount++;
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(fullyCorrectStyle);
        cell.setCellValue("the LDA/LipidBLAST columns are green, if the hit was found");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(ambigStyle);
        cell.setCellValue("the LDA/LipidBLAST column are cyan, if the hit was found additionally at a wrong RT");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(onlyMS1Style);
        cell.setCellValue("the LDA/LipidBLAST column are orange, if the hit was found only at a wrong RT");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(notFoundStyle);
        cell.setCellValue("the LDA/LipidBLAST column are red, if the hit was not reported at all or if the hit is an FP");
        
        rowCount++;
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Identification Codes:");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("2: correctly identified without any FPs");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("1: correct hit found, but FPs at wrong retention times are reported");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("0: identified only by MS1, or FP not detected");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("-1: false negative hit");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("-2: false positive hit");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("-3: false negative and false positive hit, e.g., hit at correct RT not found, but something at wrong RT reported");
      }
            
      resultWorkbook.write(out);
      out.close();
    }catch(Exception ex){
      ex.printStackTrace();
    }      
  }
  
  private CellStyle getStyleBasedOnEvidence(int evidence, boolean isFP, CellStyle fullyCorrectStyle,
      CellStyle ambigStyle, CellStyle onlyMS1Style, CellStyle notFoundStyle){
    CellStyle style=null;
    if (isFP){
      if (evidence!=LdaLBLASTCompareVO.NOT_FOUND_OR_MS1_ONLY && evidence!=LdaLBLASTCompareVO.FALSE_NEGATIVE)
        style = notFoundStyle;
    }else{
      if (evidence==LdaLBLASTCompareVO.FOUND_NO_FPS) style=fullyCorrectStyle;
      else if (evidence==LdaLBLASTCompareVO.FOUND_BUT_FPS) style=ambigStyle;
      else if (evidence==LdaLBLASTCompareVO.NOT_FOUND_OR_MS1_ONLY);//do nothing
      else style = notFoundStyle;
    }
    return style;
  }
  
  private void compareLDABLASTControlledPositiveProbes(){
    LinkedHashMap<String,LinkedHashMap<String,Boolean[]>> comparableClassesAndAdducts = new LinkedHashMap<String,LinkedHashMap<String,Boolean[]>>();
    Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>> correctAnalytes = new Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>>();
    Hashtable<String,Integer> snPositions = new Hashtable<String,Integer>();
    

    String lipidClass = "PI";
    //the first boolean tells if this adduct should be included in statistics, and the second if it the analytes should be RT-filtered
    LinkedHashMap<String,Boolean[]> adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("NH4",new Boolean[]{false,true});
    adducts.put("Na",new Boolean[]{false,true});
    adducts.put("H",new Boolean[]{false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPIStandards());
    
    lipidClass = "P-PC";
    adducts = new  LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true});
    adducts.put("Na",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPPCStandards());
    
    lipidClass = "P-PE";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true});
    adducts.put("Na",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPPEStandards());    

    lipidClass = "LPC";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true});
    adducts.put("Na",new Boolean[]{true,true});
    adducts.put("-OH",new Boolean[]{false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    correctAnalytes.put(lipidClass, getLPCStandards());

    lipidClass = "LPE";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true});
    adducts.put("Na",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    correctAnalytes.put(lipidClass, getLPEStandards());
  
    lipidClass = "PS";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,false});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPSStandards());
     
    lipidClass = "LPS";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    correctAnalytes.put(lipidClass, getLPSStandards());

    lipidClass = "PC";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true});
    adducts.put("Na",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPCStandards());

    lipidClass = "PE";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true});
    adducts.put("Na",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPEStandards());

    lipidClass = "PG";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{false,true});
    adducts.put("Na",new Boolean[]{false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPGStandards());
    
    //this has to be before DG and TG
    correctRetentionTimes(0.4d,correctAnalytes);

    lipidClass = "DG";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("Na",new Boolean[]{true,true});
    adducts.put("NH4",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 3);
    correctAnalytes.put(lipidClass, getDGStandards());

    lipidClass = "TG";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("NH4",new Boolean[]{true,true});
    adducts.put("Na",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 3);
    correctAnalytes.put(lipidClass, getTGStandards());

    lipidClass = "SM";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true});
    adducts.put("Na",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    LinkedHashMap<String,ReferenceInfoVO> smStandards = getSMStandards();
    correctRetentionTimes(0.4d,smStandards.values());
    correctAnalytes.put(lipidClass, smStandards);

    lipidClass = "Cer";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true});
    adducts.put("Na",new Boolean[]{false,true});
    adducts.put("-OH",new Boolean[]{false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    LinkedHashMap<String,ReferenceInfoVO> cerStandards = getCerStandards();
    correctRetentionTimes(0.4d,cerStandards.values());
    correctAnalytes.put(lipidClass, cerStandards);

    
    //String chromFile = "D:\\Experiment1\\LipidBLAST\\014_Ex1_Orbitrap_CID_neg_50.chrom";
    String chromFile = "D:\\Experiment1\\Orbitrap_CID\\positive\\50\\018_Ex1_Orbitrap_CID_pos_50.chrom";
    //String quantFile = "D:\\Experiment1\\massLists\\PC_PE_PG_PI_PS.xlsx";
    String quantFile = "D:\\Experiment1\\massLists\\Ex1_pos.xlsx";
    //String ldaFile = "D:\\Experiment1\\LipidBLAST\\014_Ex1_Orbitrap_CID_neg_50_PC_PE_PG_PI_PS.xlsx";
    String ldaFile = "D:\\Experiment1\\Orbitrap_CID\\positive\\50\\018_Ex1_Orbitrap_CID_pos_50_Ex1_pos.xlsx";
    //String lbFile = "D:\\Experiment1\\LipidBLAST\\output\\014_Ex1_Orbitrap_CID_neg_50.mgf_LB10.tsv";
    String lbFile = "D:\\Experiment1\\LipidBlast\\positive\\output\\018_Ex1_Orbitrap_CID_pos_50_MF450.mgf.tsv";
    //String outputFile = "D:\\Experiment1\\LipidBlast\\negative\\StandardsNegativeEvaluation_50%_LB450_comp.xlsxn";
    String outputFile = "D:\\Experiment1\\LipidBlast\\positive\\018_Ex1_Orbitrap_CID_pos_50_LB450_comp_generated.xlsx";
    compareLDABLASTControlledProbes(comparableClassesAndAdducts,snPositions,correctAnalytes,chromFile,quantFile,ldaFile,lbFile,outputFile);
  }
  
  private void compareLDABLASTControlledNegativeProbes(){
    LinkedHashMap<String,LinkedHashMap<String,Boolean[]>> comparableClassesAndAdducts = new LinkedHashMap<String,LinkedHashMap<String,Boolean[]>>();
    Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>> correctAnalytes = new Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>>();
    Hashtable<String,Integer> snPositions = new Hashtable<String,Integer>();
    

    String lipidClass = "PI";
    //the first boolean tells if this adduct should be included in statistics, and the second if it the analytes should be RT-filtered
    LinkedHashMap<String,Boolean[]> adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("-H",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPIStandards());

            
    lipidClass = "P-PC";
    adducts = new  LinkedHashMap<String,Boolean[]>();
    adducts.put("HCOO",new Boolean[]{false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPPCStandards());    

    lipidClass = "P-PE";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("-H",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPPEStandards());    

    lipidClass = "LPC";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("HCOO",new Boolean[]{false,true});
    adducts.put("-CH3",new Boolean[]{false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    correctAnalytes.put(lipidClass, getLPCStandards());

    lipidClass = "LPE";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("-H",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    correctAnalytes.put(lipidClass, getLPEStandards());
    
    lipidClass = "PS";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("-H",new Boolean[]{true,false});
    adducts.put("Na-H2",new Boolean[]{false,false});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPSStandards());

    lipidClass = "LPS";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("-H",new Boolean[]{false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    correctAnalytes.put(lipidClass, getLPSStandards());
    
    lipidClass = "PC";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("HCOO",new Boolean[]{true,true});
    adducts.put("-CH3",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPCStandards());

    lipidClass = "PE";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("-H",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPEStandards());

    lipidClass = "PG";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("-H",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    correctAnalytes.put(lipidClass, getPGStandards());

    lipidClass = "SM";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("HCOO",new Boolean[]{false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    correctAnalytes.put(lipidClass, getSMStandards());

    lipidClass = "Cer";
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("HCOO",new Boolean[]{false,true});
    adducts.put("Cl",new Boolean[]{false,true});
    adducts.put("-H",new Boolean[]{true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    correctAnalytes.put(lipidClass, getCerStandards());

    
    //String chromFile = "D:\\Experiment1\\LipidBLAST\\014_Ex1_Orbitrap_CID_neg_50.chrom";
    String chromFile = "D:\\Experiment1\\Orbitrap_CID\\negative\\50\\018_Ex1_Orbitrap_CID_neg_50.chrom";
    //String quantFile = "D:\\Experiment1\\massLists\\PC_PE_PG_PI_PS.xlsx";
    String quantFile = "D:\\Experiment1\\massLists\\Ex1_neg.xlsx";
    //String ldaFile = "D:\\Experiment1\\LipidBLAST\\014_Ex1_Orbitrap_CID_neg_50_PC_PE_PG_PI_PS.xlsx";
    String ldaFile = "D:\\Experiment1\\Orbitrap_CID\\negative\\50\\018_Ex1_Orbitrap_CID_neg_50_Ex1_neg.xlsx";
    //String lbFile = "D:\\Experiment1\\LipidBLAST\\output\\014_Ex1_Orbitrap_CID_neg_50.mgf_LB10.tsv";
    String lbFile = "D:\\Experiment1\\LipidBlast\\negative\\output\\018_Ex1_Orbitrap_CID_neg_50_MF10.mgf.tsv";
    //String outputFile = "D:\\Experiment1\\LipidBlast\\negative\\StandardsNegativeEvaluation_50%_LB450_comp.xlsxn";
    String outputFile = "D:\\Experiment1\\LipidBlast\\negative\\018_Ex1_Orbitrap_CID_neg_50_LB10_comp_generated.xlsx";
    compareLDABLASTControlledProbes(comparableClassesAndAdducts,snPositions,correctAnalytes,chromFile,quantFile,ldaFile,lbFile,outputFile);

  }
  
  private void compareLDABLASTControlledProbes(LinkedHashMap<String,LinkedHashMap<String,Boolean[]>> comparableClassesAndAdducts,
      Hashtable<String,Integer> snPositions, Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>> correctAnalytes,
      String chromFile, String quantFile, String ldaFile, String lbFile, String outputFile){
    try{
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outputFile));
      Workbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      CellStyle leftHeaderStyle = resultWorkbook.createCellStyle();
      leftHeaderStyle.cloneStyleFrom(headerStyle);
      leftHeaderStyle.setAlignment(CellStyle.ALIGN_LEFT);
      headerStyle = leftHeaderStyle;
      CellStyle centerStyle = getCenterStyle(resultWorkbook);

      CellStyle fullyCorrectStyle = getFullyCorrectStyle(resultWorkbook);
      CellStyle faCorrectStyle = getFACorrectStyle(resultWorkbook);
      CellStyle faFoundStyle = getFAFoundStyle(resultWorkbook);
      CellStyle ms1FoundStyle = getMS1FoundStyle(resultWorkbook);
      CellStyle notFoundStyle = getNotFoundStyle(resultWorkbook);
      Hashtable<String,Vector<LipidParameterSet>> resultsLDA = LDAResultReader.readResultFile(ldaFile, new Hashtable<String,Boolean>()).getIdentifications();
      LipidBLASTParser lBlastParser = new LipidBLASTParser(lbFile);
      lBlastParser.parse();
      Hashtable<String,Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>> resultsLB = lBlastParser.getResults_();
      String[] chromPaths = StringUtils.getChromFilePaths(chromFile);   
      LipidomicsAnalyzer analyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
      ////Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantVOs = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)QuantificationThread.parseQuantExcelFile(quantFile,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f).get(3);
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantVOs = null;
      Sheet summarySheet = resultWorkbook.createSheet("Summary");
      Sheet detailsSheet = resultWorkbook.createSheet("Details");

      int ms1Detectables = 0;
      int faDetectables = 0;
      int positionDetectables = 0;
      
      int ms1DetectedLDA = 0;
      int faDetectedLDAWOFPs = 0;
      int faDetectedLDA = 0;
      int posDetectedLDA = 0;

      int ms1DetectedLB = 0;
      int faDetectedLBWOFPs = 0;
      int faDetectedLB = 0;
      int posDetectedLB = 0;
      
      int rowCountDetails = 0;
      int rowCountSummary = 0;
      
      int classColumn = 0;
      int correctColumn = 1;
      int width = (12+1)*256;
      detailsSheet.setColumnWidth(correctColumn, width);
      int adductColumn = 2;
      int ldaColumn = 3;
      width = (26+1)*256;
      detailsSheet.setColumnWidth(ldaColumn, width);
      int lbColumn = 4;
      detailsSheet.setColumnWidth(lbColumn, width);
      int lbProb = 5;
      int ldaRt = 6;
/****      int lbRts = 7;
      int commentColumn = 8;*/
      int commentColumn = 7;
      
      int sumTypeColumn = 0;
      int ldaMs1Column = 1;
      int ldaMs1PercentColumn = 2;
      int lbMs1Column = 3;
      int lbMs1PercentColumn = 4;
      int ldaFaColumn = 5;
      int ldaFaPercentColumn = 6;
      int lbFaColumn = 7;
      int lbFaPercentColumn = 8;
      int ldaFaWoFPColumn = 9;
      int ldaFaWoFPPercentColumn = 10;
      int lbFaWoFPColumn = 11;
      int lbFaWoFPPercentColumn = 12;
      int ldaPosColumn = 13;
      int ldaPosPercentColumn = 14;
      int lbPosColumn = 15;
      int lbPosPercentColumn = 16;

      Row sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      Cell sumCell = sumRow.createCell(ldaMs1Column);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA identified");
      sumCell = sumRow.createCell(lbMs1Column);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LBlast identified");
      sumCell = sumRow.createCell(ldaFaColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA FAs");
      sumCell = sumRow.createCell(lbFaColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LBlast FAs");
      sumCell = sumRow.createCell(ldaFaWoFPColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA FAs no FPs");
      sumCell = sumRow.createCell(lbFaWoFPColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LBlast FAs no FPs");
      sumCell = sumRow.createCell(ldaPosColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA Pos");
      sumCell = sumRow.createCell(lbPosColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LBlast Pos");
      int typeColumn = 0;
      int ldaHitsColumn = 4;
      int ldaPercentColumn = 5;
      int lbHitsColumn = 6;
      int lbPercentColumn = 7;

      
      LinkedHashMap<String,LinkedHashMap<String,LdaLbStandardsEvidence>> speciesEvidence = new LinkedHashMap<String,LinkedHashMap<String,LdaLbStandardsEvidence>>();
      
      for (String lipidClass : comparableClassesAndAdducts.keySet()){
        LinkedHashMap<String,Boolean[]> comparableAdducts = comparableClassesAndAdducts.get(lipidClass);
        int numberSns = snPositions.get(lipidClass);
        LinkedHashMap<String,ReferenceInfoVO> ms1Analytes = correctAnalytes.get(lipidClass);
        Vector<LipidParameterSet> ldaAnalytes = resultsLDA.get(lipidClass);
        Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>> lbAnalytes = new Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>();
        if (resultsLB.containsKey(lipidClass)) lbAnalytes = resultsLB.get(lipidClass);
        
        LinkedHashMap<String,LdaLbStandardsEvidence> speciesEvidenceClass = new LinkedHashMap<String,LdaLbStandardsEvidence>();
        Row row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;
        Cell cell = row.createCell(classColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Class");
        cell = row.createCell(correctColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Correct");
        cell = row.createCell(adductColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Adduct");
        cell = row.createCell(ldaColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA");
        cell = row.createCell(lbColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LipidBLAST");
        cell = row.createCell(lbProb);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LB-Prob");
        cell = row.createCell(ldaRt);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA-RT");
/****        cell = row.createCell(lbRts);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LB-RTs");*/
        cell = row.createCell(commentColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Comments");

        int ms1DetectablesClass = 0;
        int faDetectablesClass = 0;
        int positionDetectablesClass = 0;
        
        int ms1DetectedLDAClass = 0;
        int faDetectedLDAWOFPsClass = 0;
        int faDetectedLDAClass = 0;
        int posDetectedLDAClass = 0;

        int ms1DetectedLBClass = 0;
        int faDetectedLBWOFPsClass = 0;
        int faDetectedLBClass = 0;
        int posDetectedLBClass = 0;
        
        for (String adduct : comparableAdducts.keySet()){
          
          int ms1DetectablesAdduct = 0;
          int faDetectablesAdduct = 0;
          int positionDetectablesAdduct = 0;
          
          int ms1DetectedLDAAdduct = 0;
          int faDetectedLDAWOFPsAdduct = 0;
          int faDetectedLDAAdduct = 0;
          int posDetectedLDAAdduct = 0;

          int ms1DetectedLBAdduct = 0;
          int faDetectedLBWOFPsAdduct = 0;
          int faDetectedLBAdduct = 0;
          int posDetectedLBAdduct = 0;

          
          boolean includeInAnalysis = comparableAdducts.get(adduct)[0];
          for (String analyte : ms1Analytes.keySet()){
            ReferenceInfoVO correctStructureInfo = ms1Analytes.get(analyte);
            LipidomicsMSnSet ldaAnalyte = (LipidomicsMSnSet)getLDAAnalyte(ldaAnalytes,analyte,adduct,false,true,correctStructureInfo.getCorrectRts());
            //WAS-DONE: for backward compatibility, the following was changed: getCorrectRts()[0]
            LipidBLASTIdentificationVO identVO = getLBAnalyte(lbAnalytes,analyte,adduct,correctStructureInfo.getCorrectRts()[0],comparableAdducts.get(adduct)[1]);
            String correctStructure = correctStructureInfo.getMS2Name();
            
            LdaLbStandardsEvidence ldaLbEvidence = new LdaLbStandardsEvidence(analyte);
            if (speciesEvidenceClass.containsKey(analyte)) ldaLbEvidence = speciesEvidenceClass.get(analyte);
//            for (String correctStructure : correctStructures.keySet()){
              row = detailsSheet.createRow(rowCountDetails);
              rowCountDetails++;
              cell = row.createCell(classColumn);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(lipidClass);
              cell = row.createCell(correctColumn);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(correctStructure);
              cell = row.createCell(adductColumn);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(adduct);

              Object[] lda = checkLDAEvidence(correctStructure, ldaAnalyte);
              String nameString = (String)lda[1];
              int evidence = (Integer)lda[0];
              
              Object[] lBlast = checkLBLASTEvidence(correctStructure, identVO, numberSns);
              String blastString = (String)lBlast[1];
              int blastEvidence = (Integer)lBlast[0];

              if (lipidClass.equalsIgnoreCase("TG") || (lipidClass.equalsIgnoreCase("DG")&&adduct.equalsIgnoreCase("NH4") && !nameString.equalsIgnoreCase("not reported"))){
                if (nameString.indexOf(" ")!=-1 && !lipidClass.equalsIgnoreCase("DG")) nameString = nameString.substring(0,nameString.indexOf(" "));
                if (nameString.indexOf(";")!=-1) nameString = nameString.substring(0,nameString.indexOf(";"));
                nameString = nameString.replaceAll("/", "_");
                if (evidence==FA_CORRECT) evidence = POS_CORRECT;
                if (blastEvidence==FA_CORRECT) blastEvidence = POS_CORRECT;
              }
              
              boolean msMSPresent = true;
              LipidParameterSet ldaMS1Result = null;
              if (evidence==0 && blastEvidence==0){
                ldaMS1Result = getLDAAnalyte(ldaAnalytes,analyte,adduct,true,true,correctStructureInfo.getCorrectRts());
                if (!areThereAnyMSnSpectra(analyzer,ldaMS1Result,quantVOs.get(lipidClass).get(analyte).get(adduct))){
                  evidence = -1;
                  nameString = "no MS/MS";
                  if (ldaMS1Result!=null) nameString = ldaMS1Result.getNameStringWithoutRt();      
                  blastEvidence = -1;
                  blastString = "no MS/MS";
                  msMSPresent = false;
                }               
              } 
              cell = row.createCell(ldaColumn);
              CellStyle style = getStyleBasedOnEvidence(evidence, fullyCorrectStyle, faCorrectStyle, faFoundStyle,
                  ms1FoundStyle, notFoundStyle);
              if (evidence==0 && blastEvidence>0 ){
                ldaMS1Result = getLDAAnalyte(ldaAnalytes,analyte,adduct,true,true,correctStructureInfo.getCorrectRts());
                if (ldaMS1Result!=null){
                  nameString = ldaMS1Result.getNameStringWithoutRt();
                  style=null;
                }
              }
              if (style!=null) cell.setCellStyle(style);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(nameString);

              cell = row.createCell(lbColumn);
              style = getStyleBasedOnEvidence(blastEvidence, fullyCorrectStyle, faCorrectStyle, faFoundStyle,
                  ms1FoundStyle, notFoundStyle);
              if (style!=null) cell.setCellStyle(style);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(blastString);

              if (blastEvidence>0){
                cell = row.createCell(lbProb);
                cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                cell.setCellValue((Double)lBlast[2]);
                
/****                cell = row.createCell(lbRts);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                String rtString = "";
                for (String rt : (Vector<String>)lBlast[3]) rtString+=rt+";";
                rtString = rtString.substring(0,rtString.length()-1);
                cell.setCellValue(rtString);*/
              }
              if (evidence>0 || ldaMS1Result!=null){
                cell = row.createCell(ldaRt);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                if (evidence>0)
                  cell.setCellValue(ldaAnalyte.getRt());
                else
                  cell.setCellValue(ldaMS1Result.getRt());
              }
              if (includeInAnalysis){
                if (msMSPresent && correctStructureInfo.useInEvaluation()){
                  ms1DetectablesClass++;
                  ms1DetectablesAdduct++;
                  //for debugging
//                  if ((blastEvidence>=FA_FOUND)){
//                    cell = row.createCell(commentColumn);
//                    cell.setCellType(Cell.CELL_TYPE_STRING);
//                    cell.setCellValue(((Vector<String>)lBlast[3]).toString());                
//                  }
                  if (blastEvidence==0 && identVO!=null && identVO.hasRemovedRts()){
                    cell = row.createCell(commentColumn);
                    cell.setCellType(Cell.CELL_TYPE_STRING);
                    //WAS-DONE: for backward compatibility, the following was changed: getCorrectRts()[0]
                    String startRt = Calculator.FormatNumberToString(correctStructureInfo.getCorrectRts()[0]-RT_TOL, 2);
                    String stopRt = Calculator.FormatNumberToString(correctStructureInfo.getCorrectRts()[0]+RT_TOL, 2);
                    cell.setCellValue("LB reported hits outside the RT-tolerance "+startRt+"-"+stopRt+"min: "+identVO.getOutsideRtDetections(numberSns));
                  }
                  if (evidence>=MS1_FOUND){
                    ms1DetectedLDAClass++;
                    ms1DetectedLDAAdduct++;
                    ldaLbEvidence.setMs1DetectedLDA(true);
                  }
                  if (blastEvidence>=MS1_FOUND){
                    ms1DetectedLBClass++;
                    ms1DetectedLBAdduct++;
                    ldaLbEvidence.setMs1DetectedLB(true);
                  }
                  if (numberSns>1){
                    faDetectablesClass++;
                    faDetectablesAdduct++;
                    ldaLbEvidence.setFaDetectable(true);
                    if (correctStructureInfo.isPositionAvailable() && !(lipidClass.equalsIgnoreCase("DG") && adduct.equalsIgnoreCase("NH4"))){
                      positionDetectablesClass++;
                      positionDetectablesAdduct++;
                      ldaLbEvidence.setPositionDetectable(true);
                      if (evidence>=POS_CORRECT){
                        posDetectedLDAClass++;
                        posDetectedLDAAdduct++;
                        ldaLbEvidence.setPosDetectedLDA(true);
                      }
                      if (blastEvidence>=POS_CORRECT){
                        posDetectedLBClass++;
                        posDetectedLBAdduct++;
                        ldaLbEvidence.setPosDetectedLB(true);
                      }
                    }
                    if (evidence>=FA_CORRECT){
                      faDetectedLDAWOFPsClass++;
                      faDetectedLDAWOFPsAdduct++;
                      ldaLbEvidence.setFaDetectedLDAWoFPs(true);
                    }
                    if (blastEvidence>=FA_CORRECT){
                      faDetectedLBWOFPsClass++;
                      faDetectedLBWOFPsAdduct++;
                      ldaLbEvidence.setFaDetectedLBWoFPs(true);
                    }
                    if (evidence>=FA_FOUND){
                      faDetectedLDAClass++;
                      faDetectedLDAAdduct++;
                      ldaLbEvidence.setFaDetectedLDA(true);
                    }
                    if (blastEvidence>=FA_FOUND){
                      faDetectedLBClass++;
                      faDetectedLBAdduct++;
                      ldaLbEvidence.setFaDetectedLB(true);
                    }
                  }
                  speciesEvidenceClass.put(analyte,ldaLbEvidence);
//                } else if (!correctStructures.get(correctStructure)[0]){
                }else{  
                  cell = row.createCell(commentColumn);
                  cell.setCellType(Cell.CELL_TYPE_STRING);
                  cell.setCellValue("excluded from statistics");                  
                }
              }else{
                cell = row.createCell(commentColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue("excluded from statistics");                
              }
//            }
          }
          
          if (ms1DetectablesAdduct==0) continue;
          
          sumRow = summarySheet.createRow(rowCountSummary);
          rowCountSummary++;
          sumCell = sumRow.createCell(sumTypeColumn);
          sumCell.setCellStyle(leftHeaderStyle);
          sumCell.setCellValue(lipidClass+"_"+adduct);
          
          sumCell = sumRow.createCell(ldaMs1Column);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          sumCell.setCellValue(ms1DetectedLDAAdduct+"/"+ms1DetectablesAdduct);
          sumCell = sumRow.createCell(ldaMs1PercentColumn);
          sumCell.getCellStyle().setAlignment(CellStyle.ALIGN_CENTER);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*ms1DetectedLDAAdduct))/((double)ms1DetectablesAdduct), 0)+"%");
          sumCell = sumRow.createCell(lbMs1Column);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          sumCell.setCellValue(ms1DetectedLBAdduct+"/"+ms1DetectablesAdduct);
          sumCell = sumRow.createCell(lbMs1PercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*ms1DetectedLBAdduct))/((double)ms1DetectablesAdduct), 0)+"%");
          sumCell = sumRow.createCell(ldaFaColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1) sumCell.setCellValue(faDetectedLDAAdduct+"/"+faDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(ldaFaPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLDAAdduct))/((double)faDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(lbFaColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1) sumCell.setCellValue(faDetectedLBAdduct+"/"+faDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(lbFaPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLBAdduct))/((double)faDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(ldaFaWoFPColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1) sumCell.setCellValue(faDetectedLDAWOFPsAdduct+"/"+faDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(ldaFaWoFPPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLDAWOFPsAdduct))/((double)faDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(lbFaWoFPColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1) sumCell.setCellValue(faDetectedLBWOFPsAdduct+"/"+faDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(lbFaWoFPPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLBWOFPsAdduct))/((double)faDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(ldaPosColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1 && positionDetectablesAdduct>0 && !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(posDetectedLDAAdduct+"/"+positionDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(ldaPosPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1 && positionDetectablesAdduct>0 &&  !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*posDetectedLDAAdduct))/((double)positionDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(lbPosColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1 && positionDetectablesAdduct>0 && !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(posDetectedLBAdduct+"/"+positionDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(lbPosPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (numberSns>1 && positionDetectablesAdduct>0 &&  !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*posDetectedLBAdduct))/((double)positionDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");

          
          ms1Detectables += ms1DetectablesAdduct;
          faDetectables += faDetectablesAdduct;
          positionDetectables += positionDetectablesAdduct;
          ms1DetectedLDA += ms1DetectedLDAAdduct;
          faDetectedLDAWOFPs += faDetectedLDAWOFPsAdduct;
          faDetectedLDA += faDetectedLDAAdduct;
          posDetectedLDA += posDetectedLDAAdduct;
          ms1DetectedLB += ms1DetectedLBAdduct;
          faDetectedLBWOFPs += faDetectedLBWOFPsAdduct;
          faDetectedLB += faDetectedLBAdduct;
          posDetectedLB += posDetectedLBAdduct;
          
        }
        rowCountDetails++;
        row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;
        cell = row.createCell(ldaHitsColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA");
        cell = row.createCell(lbHitsColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LipidBLAST");
                
        row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;        
        cell = row.createCell(typeColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Identified Species");
        cell = row.createCell(ldaHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        cell.setCellValue(ms1DetectedLDAClass+"/"+ms1DetectablesClass);
        cell = row.createCell(ldaPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (ms1DetectablesClass==0)
          cell.setCellValue("0%");
        else
          cell.setCellValue(Calculator.FormatNumberToString(((double)(100*ms1DetectedLDAClass))/((double)ms1DetectablesClass), 0)+"%");
        cell = row.createCell(lbHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        cell.setCellValue(ms1DetectedLBClass+"/"+ms1DetectablesClass);
        cell = row.createCell(lbPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (ms1DetectablesClass==0)
          cell.setCellValue("0%");
        else
          cell.setCellValue(Calculator.FormatNumberToString(((double)(100*ms1DetectedLBClass))/((double)ms1DetectablesClass), 0)+"%");
        
        row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;        
        cell = row.createCell(typeColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Identified Fatty acid compositions");
        cell = row.createCell(ldaHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0) cell.setCellValue(faDetectedLDAClass+"/"+faDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(ldaPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0) cell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLDAClass))/((double)faDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");
        cell = row.createCell(lbHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0) cell.setCellValue(faDetectedLBClass+"/"+faDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(lbPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0) cell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLBClass))/((double)faDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");

        row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;        
        cell = row.createCell(typeColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Identified Fatty acid compositions w.o. FPs");
        cell = row.createCell(ldaHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0) cell.setCellValue(faDetectedLDAWOFPsClass+"/"+faDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(ldaPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0) cell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLDAWOFPsClass))/((double)faDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");
        cell = row.createCell(lbHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0) cell.setCellValue(faDetectedLBWOFPsClass+"/"+faDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(lbPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0) cell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLBWOFPsClass))/((double)faDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");

        row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;        
        cell = row.createCell(typeColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Unambiguously identfied positions");
        cell = row.createCell(ldaHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) cell.setCellValue(posDetectedLDAClass+"/"+positionDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(ldaPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) cell.setCellValue(Calculator.FormatNumberToString(((double)(100*posDetectedLDAClass))/((double)positionDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");
        cell = row.createCell(lbHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1  && ms1DetectablesClass>0 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) cell.setCellValue(posDetectedLBClass+"/"+positionDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(lbPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1  && ms1DetectablesClass>0 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) cell.setCellValue(Calculator.FormatNumberToString((double)(100*posDetectedLBClass)/((double)positionDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");
        rowCountDetails++;
        rowCountDetails++;

        if (ms1DetectablesClass==0) continue;
        speciesEvidence.put(lipidClass, speciesEvidenceClass);


/****        
        sumRow = summarySheet.createRow(rowCountSummary);
        rowCountSummary++;
        sumCell = sumRow.createCell(sumTypeColumn);
        sumCell.setCellStyle(leftHeaderStyle);
        sumCell.setCellValue(lipidClass);
        
        sumCell = sumRow.createCell(ldaMs1Column);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(ms1DetectedLDAClass+"/"+ms1DetectablesClass);
        sumCell = sumRow.createCell(ldaMs1PercentColumn);
        sumCell.getCellStyle().setAlignment(CellStyle.ALIGN_CENTER);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(((100*ms1DetectedLDAClass)/ms1DetectablesClass)+"%");
        sumCell = sumRow.createCell(lbMs1Column);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(ms1DetectedLBClass+"/"+ms1DetectablesClass);
        sumCell = sumRow.createCell(lbMs1PercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(((100*ms1DetectedLBClass)/ms1DetectablesClass)+"%");
        sumCell = sumRow.createCell(ldaFaColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1) sumCell.setCellValue(faDetectedLDAClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaFaPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1) sumCell.setCellValue(((100*faDetectedLDAClass)/faDetectablesClass)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbFaColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1) sumCell.setCellValue(faDetectedLBClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbFaPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1) sumCell.setCellValue(((100*faDetectedLBClass)/faDetectablesClass)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaFaWoFPColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1) sumCell.setCellValue(faDetectedLDAWOFPsClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaFaWoFPPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1) sumCell.setCellValue(((100*faDetectedLDAWOFPsClass)/faDetectablesClass)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbFaWoFPColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1) sumCell.setCellValue(faDetectedLBWOFPsClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbFaWoFPPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1) sumCell.setCellValue(((100*faDetectedLBWOFPsClass)/faDetectablesClass)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaPosColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(posDetectedLDAClass+"/"+positionDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaPosPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1 && positionDetectablesClass>0 &&  !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(((100*posDetectedLDAClass)/positionDetectablesClass)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbPosColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(posDetectedLBClass+"/"+positionDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbPosPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (numberSns>1 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(((100*posDetectedLBClass)/positionDetectablesClass)+"%");
        else sumCell.setCellValue("NA");

        
        ms1Detectables += ms1DetectablesClass;
        faDetectables += faDetectablesClass;
        positionDetectables += positionDetectablesClass;
        ms1DetectedLDA += ms1DetectedLDAClass;
        faDetectedLDAWOFPs += faDetectedLDAWOFPsClass;
        faDetectedLDA += faDetectedLDAClass;
        posDetectedLDA += posDetectedLDAClass;
        ms1DetectedLB += ms1DetectedLBClass;
        faDetectedLBWOFPs += faDetectedLBWOFPsClass;
        faDetectedLB += faDetectedLBClass;
        posDetectedLB += posDetectedLBClass;*/
                
      }
      rowCountSummary++;
      
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("Total");
      
      sumCell = sumRow.createCell(ldaMs1Column);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(ms1DetectedLDA+"/"+ms1Detectables);
      sumCell = sumRow.createCell(ldaMs1PercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLDA)/((double)ms1Detectables), 0)+"%");
      sumCell = sumRow.createCell(lbMs1Column);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(ms1DetectedLB+"/"+ms1Detectables);
      sumCell = sumRow.createCell(lbMs1PercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLB)/((double)ms1Detectables), 0)+"%");
      sumCell = sumRow.createCell(ldaFaColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLDA+"/"+faDetectables);
      sumCell = sumRow.createCell(ldaFaPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDA)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(lbFaColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLB+"/"+faDetectables);
      sumCell = sumRow.createCell(lbFaPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLB)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(ldaFaWoFPColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLDAWOFPs+"/"+faDetectables);
      sumCell = sumRow.createCell(ldaFaWoFPPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDAWOFPs)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(lbFaWoFPColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLBWOFPs+"/"+faDetectables);
      sumCell = sumRow.createCell(lbFaWoFPPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLBWOFPs)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(ldaPosColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (positionDetectables>0) sumCell.setCellValue(posDetectedLDA+"/"+positionDetectables);
      sumCell = sumRow.createCell(ldaPosPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (positionDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*posDetectedLDA)/((double)positionDetectables), 0)+"%");
      sumCell = sumRow.createCell(lbPosColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (positionDetectables>0) sumCell.setCellValue(posDetectedLB+"/"+positionDetectables);
      sumCell = sumRow.createCell(lbPosPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (positionDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*posDetectedLB)/((double)positionDetectables), 0)+"%");
      
      rowCountSummary++;
      rowCountSummary++;
      rowCountSummary++;
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(ldaMs1Column);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA identified");
      sumCell = sumRow.createCell(lbMs1Column);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LBlast identified");
      sumCell = sumRow.createCell(ldaFaColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA FAs");
      sumCell = sumRow.createCell(lbFaColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LBlast FAs");
      sumCell = sumRow.createCell(ldaFaWoFPColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA FAs no FPs");
      sumCell = sumRow.createCell(lbFaWoFPColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LBlast FAs no FPs");
      sumCell = sumRow.createCell(ldaPosColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA Pos");
      sumCell = sumRow.createCell(lbPosColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LBlast Pos");
      
      ms1Detectables = 0;   
      faDetectables = 0;
      positionDetectables = 0;  
      ms1DetectedLDA = 0;
      faDetectedLDAWOFPs = 0;
      faDetectedLDA = 0;
      posDetectedLDA = 0;
      ms1DetectedLB = 0;
      faDetectedLBWOFPs = 0;
      faDetectedLB = 0;
      posDetectedLB = 0;
      for (String lipidClass : speciesEvidence.keySet()){
        sumRow = summarySheet.createRow(rowCountSummary);
        rowCountSummary++;
        sumCell = sumRow.createCell(sumTypeColumn);
        sumCell.setCellStyle(leftHeaderStyle);
        sumCell.setCellValue(lipidClass);
        
        LinkedHashMap<String,LdaLbStandardsEvidence> evidenceClass =  speciesEvidence.get(lipidClass);
        int ms1DetectablesClass = 0;
        int faDetectablesClass = 0;
        int positionDetectablesClass = 0;
        int ms1DetectedLDAClass = 0;
        int faDetectedLDAWOFPsClass = 0;
        int faDetectedLDAClass = 0;
        int posDetectedLDAClass = 0;
        int ms1DetectedLBClass = 0;
        int faDetectedLBWOFPsClass = 0;
        int faDetectedLBClass = 0;
        int posDetectedLBClass = 0;
        for (LdaLbStandardsEvidence ev: evidenceClass.values()){
          ms1DetectablesClass++;
          if (ev.isMs1DetectedLDA()) ms1DetectedLDAClass++;
          if (ev.isMs1DetectedLB()) ms1DetectedLBClass++;
          if (ev.isFaDetectable()){
            faDetectablesClass++;
            if (ev.isFaDetectedLDA()) faDetectedLDAClass++;
            if (ev.isFaDetectedLB()) faDetectedLBClass++;
            if (ev.isFaDetectedLDAWoFPs()) faDetectedLDAWOFPsClass++;
            if (ev.isFaDetectedLBWoFPs()) faDetectedLBWOFPsClass++;
            if (ev.isPositionDetectable()){
              positionDetectablesClass++;
              if (ev.isPosDetectedLDA()) posDetectedLDAClass++;
              if (ev.isPosDetectedLB()) posDetectedLBClass++;
            }
          }
        }
        
        sumCell = sumRow.createCell(ldaMs1Column);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(ms1DetectedLDAClass+"/"+ms1DetectablesClass);
        sumCell = sumRow.createCell(ldaMs1PercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLDAClass)/((double)ms1DetectablesClass), 0)+"%");
        sumCell = sumRow.createCell(lbMs1Column);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(ms1DetectedLBClass+"/"+ms1DetectablesClass);
        sumCell = sumRow.createCell(lbMs1PercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLBClass)/((double)ms1DetectablesClass), 0)+"%");
        sumCell = sumRow.createCell(ldaFaColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(faDetectedLDAClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaFaPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDAClass)/((double)faDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbFaColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(faDetectedLBClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbFaPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLBClass)/((double)faDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaFaWoFPColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(faDetectedLDAWOFPsClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaFaWoFPPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDAWOFPsClass)/((double)faDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbFaWoFPColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(faDetectedLBWOFPsClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbFaWoFPPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLBWOFPsClass)/((double)faDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaPosColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(posDetectedLDAClass+"/"+positionDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaPosPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0 && positionDetectablesClass>0 &&  !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*posDetectedLDAClass)/((double)positionDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbPosColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(posDetectedLBClass+"/"+positionDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(lbPosPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0 && positionDetectablesClass>0 && !lipidClass.equalsIgnoreCase("TG")) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*posDetectedLBClass)/((double)positionDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");

        ms1Detectables += ms1DetectablesClass;
        faDetectables += faDetectablesClass;
        positionDetectables += positionDetectablesClass;
        ms1DetectedLDA += ms1DetectedLDAClass;
        faDetectedLDAWOFPs += faDetectedLDAWOFPsClass;
        faDetectedLDA += faDetectedLDAClass;
        posDetectedLDA += posDetectedLDAClass;
        ms1DetectedLB += ms1DetectedLBClass;
        faDetectedLBWOFPs += faDetectedLBWOFPsClass;
        faDetectedLB += faDetectedLBClass;
        posDetectedLB += posDetectedLBClass;
      }
      rowCountSummary++;
      
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("Total");
      
      sumCell = sumRow.createCell(ldaMs1Column);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(ms1DetectedLDA+"/"+ms1Detectables);
      sumCell = sumRow.createCell(ldaMs1PercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLDA)/((double)ms1Detectables), 0)+"%");
      sumCell = sumRow.createCell(lbMs1Column);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(ms1DetectedLB+"/"+ms1Detectables);
      sumCell = sumRow.createCell(lbMs1PercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLB)/((double)ms1Detectables), 0)+"%");
      sumCell = sumRow.createCell(ldaFaColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLDA+"/"+faDetectables);
      sumCell = sumRow.createCell(ldaFaPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDA)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(lbFaColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLB+"/"+faDetectables);
      sumCell = sumRow.createCell(lbFaPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLB)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(ldaFaWoFPColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLDAWOFPs+"/"+faDetectables);
      sumCell = sumRow.createCell(ldaFaWoFPPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDAWOFPs)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(lbFaWoFPColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLBWOFPs+"/"+faDetectables);
      sumCell = sumRow.createCell(lbFaWoFPPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLBWOFPs)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(ldaPosColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (positionDetectables>0) sumCell.setCellValue(posDetectedLDA+"/"+positionDetectables);
      sumCell = sumRow.createCell(ldaPosPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (positionDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*posDetectedLDA)/((double)positionDetectables), 0)+"%");
      sumCell = sumRow.createCell(lbPosColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (positionDetectables>0) sumCell.setCellValue(posDetectedLB+"/"+positionDetectables);
      sumCell = sumRow.createCell(lbPosPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (positionDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*posDetectedLB)/((double)positionDetectables), 0)+"%");
      
      //legend
      rowCountSummary++;
      rowCountSummary++;
      rowCountSummary++;
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("Legend to the \"Details\" tab:");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(fullyCorrectStyle);
      sumCell.setCellValue("Green corresponds to a completely correct lipid species identifications");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(faCorrectStyle);
      sumCell.setCellValue("Cyan corresponds to correctly identified fatty acid (FA) position, without any false positives (FPs), and without  position assignment");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(faFoundStyle);
      sumCell.setCellValue("Dark yellow corresponds to correctly identified FA positions, but with additional FPs");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(ms1FoundStyle);
      sumCell.setCellValue("Orange corresponds to a detection of the lipid class, but without identifying the correct FAs");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(notFoundStyle);
      sumCell.setCellValue("Red corresponds to an undetected analyte (this color is not used if no MS/MS are present)");

      summarySheet.setColumnWidth(0, 4000);
      
      resultWorkbook.write(out);
      resultWorkbook.close();
      out.close();
    }catch(Exception ex){
      ex.printStackTrace();
    } 
  }
  
  private CellStyle getStyleBasedOnEvidence(int evidence, CellStyle fullyCorrectStyle, CellStyle faCorrectStyle,
      CellStyle faFoundStyle, CellStyle ms1FoundStyle, CellStyle notFoundStyle){
    CellStyle style = null;
    if (evidence==NOT_FOUND) style = notFoundStyle;
    else if (evidence==MS1_FOUND) style = ms1FoundStyle;
    else if (evidence==FA_FOUND) style = faFoundStyle;
    else if (evidence==FA_CORRECT) style = faCorrectStyle;
    else if (evidence==POS_CORRECT) style = fullyCorrectStyle;    
    return style;
  }
  
  private CellStyle getStyleBasedOnEvidenceDial(int evidence, CellStyle fullyCorrectStyle, CellStyle faCorrectStyle,
      CellStyle faFoundStyle, CellStyle ms1FoundStyle, CellStyle notFoundStyle){
    CellStyle style = null;
    if (evidence==NOT_FOUND) style = notFoundStyle;
    else if (evidence==MS1_FOUND) style = ms1FoundStyle;
    else if (evidence==FA_FOUND) style = faFoundStyle;
    else if (evidence==FA_CORRECT || evidence==POS_CORRECT) style = fullyCorrectStyle;
    return style;
  }

  
  private CellStyle getFullyCorrectStyle(Workbook wb){
    Color DARK_GREEN = new Color(0x00, 0xC0, 0x00);
    return getColorStyle(wb,DARK_GREEN);
  }
  
  private CellStyle getFACorrectStyle(Workbook wb){
    return getColorStyle(wb,Color.CYAN);
  }

  private CellStyle getFAFoundStyle(Workbook wb){
    Color VERY_DARK_YELLOW = new Color(0x80, 0x80, 0x00);
    return getColorStyle(wb,VERY_DARK_YELLOW);
  }

  private CellStyle getMS1FoundStyle(Workbook wb){
    return getColorStyle(wb,Color.ORANGE);
  }
  
  private CellStyle getNotFoundStyle(Workbook wb){
    return getColorStyle(wb,Color.RED);
  }
  
  private CellStyle getColorStyle(Workbook wb, Color awtColor){
    CellStyle style = wb.createCellStyle();
    XSSFColor color = new XSSFColor(awtColor);
    Font font = wb.createFont();
    ((XSSFFont)font).setColor(color);
    style.setFont(font);
    return style;
  }
  
  private Object[] checkLBLASTEvidence(String correctStructure, LipidBLASTIdentificationVO identVO, int amountSns){
    int evidence = 0;
    String name = "not reported";
    double maxProb = 0d;
    Vector<String> allRts = new Vector<String>();
    if (identVO!=null && identVO.hasHitsWithinRt()){
      evidence = 1;
      boolean foundCorrectFAs = false;
      boolean foundWrongFAs = false;
      boolean foundCorrectPosition = false;
      boolean foundWrongPosition = false;
      Vector<String> idents = identVO.getIdentifications(amountSns);
      String proposedDisplayString = "";
      for (String ident : idents){
        proposedDisplayString += ident+" ";
        boolean[] correctness = checkForFACorrectness(correctStructure,ident);
//        if (!foundCorrectFAs || (foundCorrectFAs && !foundCorrectPosition && correctness[1])){
//          maxProb = 0d;
//          allRts = new Vector<String>();
//        }
//        if ((correctness[0] && !foundCorrectPosition) || correctness[1]){
//          double prob = identVO.getHighestProbability(ident,amountSns);
//          Vector<String> rts = identVO.getRetentionTimes(ident,amountSns);
//          if (prob>maxProb) maxProb = prob;
//          allRts.addAll(rts);
//        } 
        if (!foundCorrectFAs || correctness[0]){
          maxProb = 0d;
          allRts = new Vector<String>();
        }
        if (!foundCorrectFAs || correctness[0]){
          double prob = identVO.getHighestProbability(ident,amountSns);
          Vector<String> rts = identVO.getRetentionTimes(ident,amountSns);
          if (prob>maxProb) maxProb = prob;
          allRts.addAll(rts);
        } 
        
        if (correctness[0]) foundCorrectFAs = true;
        else foundWrongFAs = true;
        if (correctness[1]) foundCorrectPosition = true;
        else foundWrongPosition = true;
      }
      if (proposedDisplayString.length()>0){
        proposedDisplayString.substring(0, proposedDisplayString.length()-1);
        name = proposedDisplayString;
      }
      if (foundCorrectFAs){
        evidence = 2;
        if (!foundWrongFAs){
          evidence = 3;
          if (foundCorrectPosition && !foundWrongPosition)
            evidence = 4;
        }
      }
    }
    Object[] result = new Object[4];
    result[0] = evidence;
    result[1] = name;
    result[2] = maxProb;
    result[3] = allRts;
    return result;
  }
  
  private Object[] checkDialEvidence(String lClass, String ms1Name, String correctStructure, Vector<MSDialEntry> entries, int amountSns){
    int evidence = 0;
    String name = "not reported";
    double maxProb = 0d;
    double maxProbMS1 = 0;
    Vector<String> allRts = new Vector<String>();
    boolean foundSumSpecies = false;
    boolean foundSomethingElse = false;
    boolean foundCorrectFAs = false;
    boolean foundWrongFAs = false;
    String proposedDisplayString = "";

    for (MSDialEntry entry : entries) {
      if (lClass.equalsIgnoreCase(entry.getLdaClassName()) && ms1Name.equalsIgnoreCase(entry.getLdaMs1Name())) {
        foundSumSpecies = true;
        proposedDisplayString += entry.getDialClassName()+" "+(entry.getDialMs2Name()==null ? entry.getDialMs1Name() : entry.getDialMs2Name())+" ";
        if (entry.getLdaMs2Name()!=null) {
          boolean[] correctness = checkForFACorrectness(correctStructure,entry.getLdaMs2Name());
          if (!foundCorrectFAs || correctness[0]){
            maxProb = 0d;
            allRts = new Vector<String>();
          }
          if (!foundCorrectFAs || correctness[0]){
            double prob = (double)entry.getScore();
            if (prob>maxProb) maxProb = prob;
          }
          allRts.add(Calculator.FormatNumberToString(entry.getRt(),2d));
          if (correctness[0]) foundCorrectFAs = true;
          else foundWrongFAs = true;
        } else if (entry.getLdaMs1Name()!=null) {
          double prob = (double)entry.getScore();
          if (prob>maxProbMS1) maxProbMS1 = prob;
        }
      }else {
        proposedDisplayString += entry.getDialClassName()+" "+(entry.getDialMs2Name()==null ? entry.getDialMs1Name() : entry.getDialMs2Name())+" ";
        foundSomethingElse = true;
//        System.out.println(lClass+" "+correctStructure);
//        System.out.println("Wrong: "+entry.getName());
      }
      if (proposedDisplayString.length()>0){
        proposedDisplayString.substring(0, proposedDisplayString.length()-1);
        name = proposedDisplayString;
      }
    }
    if (foundSumSpecies) {
      evidence = 1;
      if (foundCorrectFAs){
        evidence = 2;
        if (!foundWrongFAs){
          evidence = 3;
        }
      }
    }
    Object[] result = new Object[4];
    result[0] = evidence;
    result[1] = name;
    result[2] = maxProb;
    if (maxProb<=0)
      result[2] = maxProbMS1;
    result[3] = allRts;
    return result;
  }


  
  private Object[] checkLDAEvidence(String correct, LipidomicsMSnSet lda) throws LipidCombinameEncodingException{
    int evidence = 0;
    String name = "not reported";
    if (lda!=null){
      evidence = 1;
      name = lda.getNameStringWithoutRt();
      String proposedDisplayString = "";
      boolean foundCorrectFAs = false;
      boolean foundWrongFAs = false;
      boolean foundCorrectPosition = false;
      boolean foundWrongPosition = false;
      for (Object names : lda.getMSnIdentificationNames()){
        if (names instanceof Vector){
          String nameString = "";
          for (String faIds : (Vector<String>)names){
            nameString += faIds+";";
            boolean[] correctness = checkForFACorrectness(correct,faIds);
            if (correctness[0]) foundCorrectFAs = true;
            else foundWrongFAs = true;
            if (correctness[1]) foundCorrectPosition = true;
            else foundWrongPosition = true;
          }
          nameString = nameString.substring(0,nameString.length()-1);
          proposedDisplayString += nameString+" ";
        }else{
          String faIds = (String)names;
          proposedDisplayString += faIds+" ";
          boolean[] correctness = checkForFACorrectness(correct,faIds);
          if (correctness[0]) foundCorrectFAs = true;
          else foundWrongFAs = true;
          if (correctness[1]) foundCorrectPosition = true;
          else foundWrongPosition = true;
        }
      }
      if (proposedDisplayString.length()>0){
        proposedDisplayString.substring(0, proposedDisplayString.length()-1);
        name = proposedDisplayString;
      }
      if (foundCorrectFAs){
        evidence = 2;
        if (!foundWrongFAs){
          evidence = 3;
          if (foundCorrectPosition && !foundWrongPosition)
            evidence = 4;
        }
      }
    }
    Object[] result = new Object[2];
    result[0] = evidence;
    result[1] = name;
    return result;
  }
  
  private boolean[] checkForFACorrectness(String correct, String faIds){
    boolean correctFA = false;
    boolean correctPos = false;
    if (correct.equalsIgnoreCase(faIds)){
      correctFA = true;
      correctPos = true;
    } else{
      Hashtable<String,Integer> amountsCorrect = getAmountFAs(correct);
      Hashtable<String,Integer> amountsCurrent = getAmountFAs(faIds);
      boolean allCorrect = true;
      for (String fa: amountsCurrent.keySet()){
        if (!amountsCorrect.containsKey(fa) || amountsCorrect.get(fa)!=amountsCurrent.get(fa)){
          allCorrect = false;
          break;
        }
      }
      if (allCorrect) correctFA = true;
    }
    boolean[] result = new boolean[2];
    result[0] = correctFA;
    result[1] = correctPos;
    return result;
  }
  
  private Hashtable<String,Integer> getAmountFAs(String faString){
    String[] fas = faString.replaceAll("/", "_").split("_");
    Hashtable<String,Integer> amounts = new Hashtable<String,Integer>();
    for (int i=0;i!=fas.length;i++){
      int amount = 0;
      if (amounts.containsKey(fas[i])) amount = amounts.get(fas[i]);
      amount++;
      amounts.put(fas[i],amount);
    }
    return amounts;
  }
  
  private LipidParameterSet getLDAAnalyte(Vector<LipidParameterSet> ldaAnalytes, String analyte, String adduct, boolean ignoreMSn, boolean ignoreRt, double[] rts){
    LipidParameterSet result = null;
    Vector<LipidParameterSet> sets = getLDAAnalytes(ldaAnalytes, analyte, adduct, ignoreMSn, ignoreRt, rts);
    for (LipidParameterSet set : sets){
      if (result==null || result.Area<set.Area) result = set;
    }
    return result;
  }
  
  private  Vector<LipidParameterSet> getLDAAnalytes(Vector<LipidParameterSet> ldaAnalytes, String analyte, String adduct, boolean ignoreMSn){
    return getLDAAnalytes(ldaAnalytes, analyte, adduct, ignoreMSn, true, new double[0]);
  }
  
  
  private  Vector<LipidParameterSet> getLDAAnalytes(Vector<LipidParameterSet> ldaAnalytes, String analyte, String adduct, boolean ignoreMSn, boolean ignoreRt, double[] rts){
    Vector<LipidParameterSet> results = new Vector<LipidParameterSet>();
    for (LipidParameterSet set : ldaAnalytes){
      if (!set.getNameStringWithoutRt().equalsIgnoreCase(analyte)) continue;
      if (!set.getModificationName().equalsIgnoreCase(adduct)) continue;
      if (!ignoreMSn && !(set instanceof LipidomicsMSnSet)) continue;
      boolean rtOK = false;
      double ldaRt = Double.parseDouble(set.getRt());
      for (double rt : rts) {
        if ((rt-MSDIAL_RT_TOL)<ldaRt && ldaRt<(rt+MSDIAL_RT_TOL))
          rtOK = true;
      }
      if (rtOK||ignoreRt)
        results.add(set);
    }
    return results;
  }
  
  
  private LipidBLASTIdentificationVO getLBAnalyte(Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>> analytes, String analyte, String adduct, double rt, boolean doRtFiltering){
    LipidBLASTIdentificationVO result = null;
    if (analytes.containsKey(analyte) && analytes.get(analyte).containsKey(adduct)){
      result = analytes.get(analyte).get(adduct);
      if (doRtFiltering) result.filterAwayWrongRts(rt,RT_TOL);
    }
    return result;
  }
  
  private Vector<Vector<MSDialEntry>> getDialAnalytes(Vector<MSDialEntry> dialHits, double mz, double[] rts, boolean doRtFiltering){
    Vector<MSDialEntry> result = new Vector<MSDialEntry>();
    Vector<MSDialEntry> outsideRt = new Vector<MSDialEntry>();
    for (MSDialEntry entry : dialHits) {
      if ((mz-MSDIAL_MZ_TOL)<entry.getMz() && entry.getMz()<(mz+MSDIAL_MZ_TOL)){
        boolean rtOK = false;
        if (doRtFiltering) {
          for (int i=0; i!=rts.length; i++) {
            float rt = (float)rts[i];
            if (((rt-MSDIAL_RT_TOL)<entry.getRt() && entry.getRt()<(rt+MSDIAL_RT_TOL)))
              rtOK = true;
          }
        }else
          rtOK = true;
        if (rtOK)
          result.add(entry);
        else
          outsideRt.add(entry);
      }
    }
    Vector<Vector<MSDialEntry>> both = new Vector<Vector<MSDialEntry>>();
    both.add(result);
    both.add(outsideRt);
    return both;
  }
  
  private LinkedHashMap<String,ReferenceInfoVO> getLPCStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("13:0", new ReferenceInfoVO("13:0",new double[]{1.7d},true,false));
    standards.put("14:0", new ReferenceInfoVO("14:0",new double[]{2.3d},true,false));
    standards.put("16:0", new ReferenceInfoVO("16:0",new double[]{4.1d},true,false));
    standards.put("18:0", new ReferenceInfoVO("18:0",new double[]{6.8d},true,false));
    standards.put("18:1", new ReferenceInfoVO("18:1",new double[]{4.8d},true,false));    
    return standards;
  }
  
  private LinkedHashMap<String,ReferenceInfoVO> getLPCExp2Standards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("13:0", new ReferenceInfoVO("13:0",new double[]{1.7d},true,false));
    standards.put("15:0", new ReferenceInfoVO("15:0",new double[]{2.3d},true,false));
    return standards;
  }

  
  private LinkedHashMap<String,ReferenceInfoVO> getLPEStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("14:0", new ReferenceInfoVO("14:0",new double[]{2.4d},true,false));
    standards.put("18:0", new ReferenceInfoVO("18:0",new double[]{7.1d},true,false));
    standards.put("18:1", new ReferenceInfoVO("18:1",new double[]{5.0d},true,false));    
    return standards;
  }
  
  private LinkedHashMap<String,ReferenceInfoVO> getLPEExp2Standards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("16:0", new ReferenceInfoVO("16:0",new double[]{2.4d},true,false));
    standards.put("18:0", new ReferenceInfoVO("18:0",new double[]{7.1d},true,false));
    return standards;
  }

  
  private LinkedHashMap<String,ReferenceInfoVO> getPSStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("32:0", new ReferenceInfoVO("16:0/16:0",new double[]{24.8d},true,false));
    standards.put("34:0", new ReferenceInfoVO("17:0/17:0",new double[]{27.2d},true,false));
    standards.put("34:1", new ReferenceInfoVO("16:0/18:1",new double[]{24.8d},true,true));
    standards.put("34:2", new ReferenceInfoVO("16:0/18:2",new double[]{23.0d},true,true));
    standards.put("36:0", new ReferenceInfoVO("18:0/18:0",new double[]{31.0d},true,false));
    standards.put("36:1", new ReferenceInfoVO("18:0/18:1",new double[]{29.0d},true,true));
    standards.put("36:2", new ReferenceInfoVO("18:0/18:2",new double[]{26.0d},true,true));
    return standards;
  }
  
  private LinkedHashMap<String,ReferenceInfoVO> getLPSStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("16:0", new ReferenceInfoVO("16:0",new double[]{3.6d},true,false));
    standards.put("17:1", new ReferenceInfoVO("17:1",new double[]{3.2d},true,false));
    standards.put("18:0", new ReferenceInfoVO("18:0",new double[]{6.5d},true,false));
    standards.put("18:1", new ReferenceInfoVO("18:1",new double[]{4.4d},true,false));    
    return standards;
  }

  private LinkedHashMap<String,ReferenceInfoVO> getPCStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("32:0", new ReferenceInfoVO("14:0/18:0",new double[]{25.3d},true,true));
    standards.put("32:1", new ReferenceInfoVO("18:1/14:0",new double[]{23.3d},true,true));
    standards.put("34:0", new ReferenceInfoVO("16:0/18:0",new double[]{27.6d},true,true));
    standards.put("34:1", new ReferenceInfoVO("16:0/18:1",new double[]{25.9d},true,true));
    standards.put("36:0", new ReferenceInfoVO("18:0/18:0",new double[]{29.9d},true,false));
    standards.put("36:1", new ReferenceInfoVO("18:0/18:1",new double[]{28.3d},true,true));
    standards.put("36:2", new ReferenceInfoVO("18:0/18:2",new double[]{26.6d},true,true));
    standards.put("40:0", new ReferenceInfoVO("20:0/20:0",new double[]{33.7d},true,false));
    standards.put("48:2", new ReferenceInfoVO("24:1/24:1",new double[]{36.8d},true,false));
    return standards;
  }

  private LinkedHashMap<String,ReferenceInfoVO> getPCExp2Standards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("31:1", new ReferenceInfoVO("17:0/14:1",new double[]{23.3d},true,true));
    standards.put("40:6", new ReferenceInfoVO("18:0/22:6",new double[]{33.7d},true,false));
    return standards;
  }

  private LinkedHashMap<String,ReferenceInfoVO> getPEStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("32:0", new ReferenceInfoVO("16:0/16:0",new double[]{25.7d},true,false));
    standards.put("34:0", new ReferenceInfoVO("17:0/17:0",new double[]{28.1d},true,false));
    standards.put("34:1", new ReferenceInfoVO("16:0/18:1",new double[]{26.3d},true,true));
    standards.put("36:0", new ReferenceInfoVO("18:0/18:0",new double[]{30.3d},true,false));
    standards.put("36:1", new ReferenceInfoVO("18:0/18:1",new double[]{28.7d},true,true));
    standards.put("36:2", new ReferenceInfoVO("18:0/18:2",new double[]{27.0d},true,true));
    standards.put("36:4", new ReferenceInfoVO("16:0/20:4",new double[]{24.5d},true,true));
    return standards;
  }
  
  private LinkedHashMap<String,ReferenceInfoVO> getPEExp2Standards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("34:1", new ReferenceInfoVO("16:0/18:1",new double[]{26.3d},true,true));
    standards.put("43:6", new ReferenceInfoVO("21:0/22:6",new double[]{24.5d},true,true));
    return standards;
  }

  
  private LinkedHashMap<String,ReferenceInfoVO> getPIStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("25:0", new ReferenceInfoVO("12:0/13:0",new double[]{12.6d},true,true));
    standards.put("31:1", new ReferenceInfoVO("17:0/14:1",new double[]{19.7d},true,true));
    standards.put("37:4", new ReferenceInfoVO("17:0/20:4",new double[]{23.1d},true,true));
    standards.put("43:6", new ReferenceInfoVO("21:0/22:6",new double[]{27.4d},true,true));
    return standards;
  }

  private LinkedHashMap<String,ReferenceInfoVO> getPPCStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("36:1", new ReferenceInfoVO("P-18:0/18:1",new double[]{29.3d},true,false));
    standards.put("38:4", new ReferenceInfoVO("P-18:0/20:4",new double[]{27.6d},true,false));
    standards.put("40:6", new ReferenceInfoVO("P-18:0/22:6",new double[]{27.6d},true,false));    
    return standards;
  }

  private LinkedHashMap<String,ReferenceInfoVO> getPPEStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("36:1", new ReferenceInfoVO("P-18:0/18:1",new double[]{29.8d},true,false));
    standards.put("38:4", new ReferenceInfoVO("P-18:0/20:4",new double[]{28.1d},true,false));
    standards.put("40:6", new ReferenceInfoVO("P-18:0/22:6",new double[]{27.5d},true,false));    
    return standards;
  }

  private LinkedHashMap<String,ReferenceInfoVO> getPGStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("32:0", new ReferenceInfoVO("16:0/16:0",new double[]{23.8d},true,false));
    standards.put("34:0", new ReferenceInfoVO("17:0/17:0",new double[]{26.2d},true,false));
    standards.put("34:1", new ReferenceInfoVO("16:0/18:1",new double[]{24.4d},true,true));
    standards.put("36:0", new ReferenceInfoVO("18:0/18:0",new double[]{28.4d},true,false));
    standards.put("36:1", new ReferenceInfoVO("18:0/18:1",new double[]{26.8d},true,true));
    standards.put("36:2", new ReferenceInfoVO("18:0/18:2",new double[]{25.1d},true,true));
    return standards;
  }
  
  private LinkedHashMap<String,ReferenceInfoVO> getSMStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("18:0", new ReferenceInfoVO("18:0",new double[]{24.6d},true,false));
    standards.put("18:1", new ReferenceInfoVO("18:1",new double[]{22.7d},true,false));    
    standards.put("24:0", new ReferenceInfoVO("24:0",new double[]{31.8d},true,false));
    standards.put("24:1", new ReferenceInfoVO("24:1",new double[]{29.6d},true,false));
    return standards;
  }

  private LinkedHashMap<String,ReferenceInfoVO> getCerStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("16:0", new ReferenceInfoVO("16:0",new double[]{25.5d},true,false));
    standards.put("17:0", new ReferenceInfoVO("17:0",new double[]{26.8d},true,false));
    standards.put("18:0", new ReferenceInfoVO("18:0",new double[]{28.0d},true,false));
    standards.put("18:1", new ReferenceInfoVO("18:1",new double[]{26.3d},true,false));    
    standards.put("20:0", new ReferenceInfoVO("20:0",new double[]{30.4d},true,false));
    return standards;
  }

  
  private LinkedHashMap<String,ReferenceInfoVO> getDGStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("24:0", new ReferenceInfoVO("12:0/12:0/-",new double[]{20.1d},true,true));
    standards.put("32:0", new ReferenceInfoVO("16:0/16:0/-",new double[]{30.5d},true,true));
    standards.put("34:0", new ReferenceInfoVO("18:0/16:0/-",new double[]{32.6d},true,true));
    standards.put("34:1", new ReferenceInfoVO("16:0/-/18:1",new double[]{30.9d},true,true));
    standards.put("36:0", new ReferenceInfoVO("18:0/18:0/-",new double[]{34.4d},true,true));
    standards.put("36:2", new ReferenceInfoVO("18:1/-/18:1",new double[]{31.9d},true,true));
    standards.put("36:4", new ReferenceInfoVO("18:2/-/18:2",new double[]{28.2d},true,true));
    standards.put("38:0", new ReferenceInfoVO("20:0/18:0/-",new double[]{35.9d},true,true));
    return standards;
  }
  
  private LinkedHashMap<String,ReferenceInfoVO> getTGStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d44:1", new ReferenceInfoVO("14:0/16:1/14:0",new double[]{38.4d},false,false));
    standards.put("48:0", new ReferenceInfoVO("16:0/16:0/16:0",new double[]{41.6d},true,false));
    standards.put("d48:1", new ReferenceInfoVO("15:0/18:1/15:0",new double[]{40.7d},false,false));
    standards.put("50:0", new ReferenceInfoVO("16:0/16:0/18:0",new double[]{42.7d},true,false));
    standards.put("50:1", new ReferenceInfoVO("16:0/16:0/18:1",new double[]{41.8d},true,false));
    standards.put("d51:1", new ReferenceInfoVO("17:0/17:1/17:0",new double[]{42.2d},false,false));
    standards.put("d58:7", new ReferenceInfoVO("20:2/18:3/20:2",new double[]{40.4d},false,false));
    standards.put("d58:10", new ReferenceInfoVO("20:4/18:2/20:4",new double[]{38.8d},false,false));
    standards.put("d60:1", new ReferenceInfoVO("20:0/20:1/20:0",new double[]{38.8d},false,false));
    standards.put("d62:16", new ReferenceInfoVO("20:5/22:6/20:5",new double[]{36.1d},false,false));
    return standards;
  }
  
  private boolean areThereAnyMSnSpectra(LipidomicsAnalyzer analyzer,LipidParameterSet lda,QuantVO quantVO) throws CgException{
    boolean spectraThere = false;
    if (lda!=null){
      Hashtable<Integer,Boolean> msLevels =  analyzer.prepareMSnSpectraCache((float)(lda.Mz[0]-LipidomicsConstants.getMs2PrecursorTolerance()), (float)(lda.Mz[0]+LipidomicsConstants.getMs2PrecursorTolerance()),100);
      for (CgProbe probe : lda.getIsotopicProbes().get(0)){
        if (analyzer.areMSnSpectraInThisRegion(probe.LowerValley,probe.UpperValley,2)){
          spectraThere = true;
          break;
        }
      }
    }else{
      Hashtable<Integer,Boolean> msLevels =  analyzer.prepareMSnSpectraCache((float)(quantVO.getAnalyteMass()-LipidomicsConstants.getMs2PrecursorTolerance()), (float)(quantVO.getAnalyteMass()+LipidomicsConstants.getMs2PrecursorTolerance()),100);
      if (analyzer.areMSnSpectraInThisRegion(0f,Float.MAX_VALUE,2)){
        spectraThere = true;
      }
    }
    return spectraThere;
  }

  private boolean areThereAnyMSnSpectra(LipidomicsAnalyzer analyzer, float mz, LipidClassInfoVO info, LinkedHashMap<String,ReferenceInfoVO> species) throws CgException{
    boolean spectraThere = false;
    Hashtable<Integer,Boolean> msLevels =  analyzer.prepareMSnSpectraCache(mz-LipidomicsConstants.getMs2PrecursorTolerance(), mz+LipidomicsConstants.getMs2PrecursorTolerance(),100);
    if (info.checkRt() && msLevels!=null && msLevels.containsKey(2) && msLevels.get(2)){
      for (ReferenceInfoVO ref : species.values()){
        for (int i=0; i!=ref.getCorrectRts().length; i++) {
          float startRt  = (float)((ref.getCorrectRts()[i]-info.getRtTolerance())*60d);
          float stopRt  = (float)((ref.getCorrectRts()[i]+info.getRtTolerance())*60d);
          if (analyzer.areMSnSpectraInThisRegion(startRt,stopRt,2)){
            spectraThere = true;
            break;
          }
        }
        if (spectraThere)
          break;
      }
      
    } else {
      if (analyzer.areMSnSpectraInThisRegion(0f,Float.MAX_VALUE,2)) spectraThere = true;
    }
    return spectraThere;
  }
  
  private LdaLBLASTCompareVO getLdaLBlastSpeciesComparision(String ms1Analyte, LipidClassInfoVO info, LinkedHashMap<String,ReferenceInfoVO> species,
      Vector<LipidParameterSet> ldaAnalytes, LipidBLASTIdentificationVO ident,LipidomicsAnalyzer analyzer,QuantVO quant) throws Exception{    
    //split LDA identifications in MS1 identifications and ones where spectra were found
    //this vector contains identifications where no ms2 is present
    Vector<LipidParameterSet> ldaMS1Only = new Vector<LipidParameterSet>();
    //this vector contains all ms2 identifications, includding the ones where no structural information can be obtained
    Vector<LipidParameterSet> ldaMs2SpectraIdentified = new Vector<LipidParameterSet>();
    //this hashtable contains identifications where the fatty acids were determined
    Hashtable<String,Vector<LipidomicsMSnSet>> ldaMs2Evidence = new Hashtable<String,Vector<LipidomicsMSnSet>>();
    Hashtable<String,String> ldaOKAccordingToReference = new Hashtable<String,String>();
    Hashtable<String,String> lbOKAccordingToReference = new Hashtable<String,String>();

    boolean isFP = false;
    for (LipidParameterSet set : ldaAnalytes){
      if (set instanceof LipidomicsMSnSet){
        LipidomicsMSnSet ldaAnalyte = (LipidomicsMSnSet)set;
        ldaMs2SpectraIdentified.add(set);
        boolean ms1Only = ldaAnalyte.getStatus()==2 ? true : false;
        if (!ms1Only) {
          for (Object msnNames : ldaAnalyte.getMSnIdentificationNames()){
            Vector<LipidomicsMSnSet> sameMs2Ident = new Vector<LipidomicsMSnSet>();
            String nameString = null;
            if (msnNames instanceof Vector){
              nameString = ((Vector<String>)msnNames).get(0);
            }else{
              nameString = ((String)msnNames);
            }
            nameString = nameString.replaceAll("/", "_");
            ////if (quant.getAnalyteClass().equalsIgnoreCase("DG")&&set.getModificationName().equalsIgnoreCase("NH4"))nameString+="_-";
            for (String key : ldaMs2Evidence.keySet()){
              if (StaticUtils.isAPermutedVersion(key, nameString,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
                nameString = key;
                sameMs2Ident = ldaMs2Evidence.get(key);
              }
            }
            sameMs2Ident.add(ldaAnalyte);
            ldaMs2Evidence.put(nameString, sameMs2Ident);
          }
        }
      } else ldaMS1Only.add(set);
    }
    
    String correctMS1Name = new String(ms1Analyte);
    int ldaMS1Evidence = 0;
    String ldaMS1Name = "not reported";
    int lbMS1Evidence = 0;
    String lbMS1Name = "not reported";
    String rtIdeal = "";
    String rtLDA = "";
    String rtLB = "";
    LinkedHashMap<String,String> rtLBs = new LinkedHashMap<String,String>();
    double lbProb = 0d;
    boolean noMS2 = false;
    boolean ms1Only = false;
    Vector<LdaLBLASTCompareVO> ms2Compare = new Vector<LdaLBLASTCompareVO>();
    //now check if LDA or LB identified the analyte at the correct RT
    if (species!=null && species.size()>0){
      boolean ms1AnalyteDetected = false;
      //WAS-DONE: for backward compatibility, the following was changed: getCorrectRts()[0]
      rtIdeal =  Calculator.FormatNumberToString(species.values().iterator().next().getCorrectRts()[0],1);
      // this is the MS1 check of the LDA
      ldaMS1Evidence = -1;
      for (LipidParameterSet set : ldaMs2SpectraIdentified){
        boolean isAtCorrectRt = isHitAtCorrectRt(Double.parseDouble(set.getRt()), info, species);
        if (isAtCorrectRt){
          ms1AnalyteDetected = true;
          if (ldaMS1Evidence==-3 || ldaMS1Evidence==-2) ldaMS1Evidence=1;
          else if (ldaMS1Evidence==-1 || ldaMS1Evidence==0) ldaMS1Evidence=2;
        } else {
          if (ldaMS1Evidence==-2||ldaMS1Evidence==-1) ldaMS1Evidence=-3;
          else if (ldaMS1Evidence==0) ldaMS1Evidence=-2;
          else if (ldaMS1Evidence==2) ldaMS1Evidence=1;          
        }
        rtLDA += set.getRt()+" ";
      }
      if (rtLDA.length()>0) rtLDA = rtLDA.substring(0,rtLDA.length()-1);      
      //this is the ms1 check of LipidBlast
      lbMS1Evidence = -1;
      if (ident!=null){
        for (LipidBLASTDetectionVO detect : ident.getDetections()){
          rtLBs.put(detect.getRetentionTime(), detect.getRetentionTime());
          boolean isAtCorrectRt = isHitAtCorrectRt(Double.parseDouble(detect.getRetentionTime()), info, species);
          if (isAtCorrectRt){
            ms1AnalyteDetected = true;
            if (detect.getProbability()>lbProb) lbProb = detect.getProbability();
            if (lbMS1Evidence==-3 || lbMS1Evidence==-2) lbMS1Evidence=1;
            else if (lbMS1Evidence==-1 || lbMS1Evidence==0) lbMS1Evidence=2;            
          } else {
            if (lbMS1Evidence==-2||lbMS1Evidence==-1) lbMS1Evidence=-3;
            else if (lbMS1Evidence==0) lbMS1Evidence=-2;
            else if (lbMS1Evidence==2) lbMS1Evidence=1;          
          }        
        }
      }
      for (String lbRt : rtLBs.keySet()) rtLB += lbRt+" ";
      if (rtLB.length()>0) rtLB = rtLB.substring(0,rtLB.length()-1);
      if ((ldaMS1Evidence==-1||ldaMS1Evidence==-3)  && (lbMS1Evidence==-1||lbMS1Evidence==-3)){
        if (!areThereAnyMSnSpectra(analyzer,(float)quant.getAnalyteMass(),info,species)){
          if (ldaMS1Evidence==-3) ldaMS1Evidence=-2;
          if (lbMS1Evidence==-3) lbMS1Evidence=-2;
          if (ldaMS1Evidence==-1) ldaMS1Evidence=0;
          if (lbMS1Evidence==-1) lbMS1Evidence=0;
          noMS2 = true;
        }else{
          ms1AnalyteDetected = true;
        }
      }
      
      //this sets the displayed names for the MS1 identifications for LDA and LipidBlast
      if (ldaMS1Evidence>0 || ldaMS1Evidence<-1) ldaMS1Name = ms1Analyte;
      else{
        boolean foundMS1Only = false;
        String rt = "";
        for (LipidParameterSet set : ldaMS1Only){
          if (isHitAtCorrectRt(Double.parseDouble(set.getRt()), info, species)){
            rt += set.getRt()+" ";
            foundMS1Only = true;
          }
        }
        if (foundMS1Only){
          ldaMS1Name = ms1Analyte;
          ms1Only = true;
          rt = rt.substring(0,rt.length()-1);
          rtLDA = rt;
          if (lbMS1Evidence<1 && ldaMS1Evidence==-1){
            ms1AnalyteDetected = true;
            ldaMS1Evidence=0;
          }
        }
      }
      if (lbMS1Evidence>0 || lbMS1Evidence<-1) lbMS1Name = ms1Analyte;
      if (!ms1AnalyteDetected && !noMS2){
        correctMS1Name += "_FP";
        isFP = true;
      }
      
      //here the check for the ms2 evidence starts
      if (info.getSns()>1 && ms1AnalyteDetected && (ldaMS1Evidence>0||lbMS1Evidence>0)){
        for (String key : species.keySet()){
          ReferenceInfoVO ref = species.get(key);
          String keyWOSlash = key.replaceAll("/", "_");
          boolean ms2AnalyteDetected = false;
          String ldaMS2Name = "not reported";
          //WAS-DONE: for backward compatibility, the following was changed: getCorrectRts()[0]
          String rtMS2Ideal = Calculator.FormatNumberToString(ref.getCorrectRts()[0],1);
          String rtMS2LDA = "";
          String rtMS2LB = "";
          // here starts the LDA check of the MS2 evidence
          int ldaMS2Evidence = -1;
          String nameString = null;
          for (String ldaName : ldaMs2Evidence.keySet()){
            Vector<LipidomicsMSnSet> ldaHits = ldaMs2Evidence.get(ldaName);
            if (ldaHits.size()==0) continue;
            String toCompare = new String(ldaName);
            if (quant.getAnalyteClass().equalsIgnoreCase("DG")&&ldaHits.get(0).getModificationName().equalsIgnoreCase("NH4")){
              toCompare += "_-";
            }
                                                                                                                             
            if (!StaticUtils.isAPermutedVersion(keyWOSlash, toCompare,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) continue;
            ldaOKAccordingToReference.put(ldaName, ldaName);
            boolean nameDidNotChange = false;
            for (LipidomicsMSnSet set : ldaHits){              
              // display name for LDA
              for (Object msnNames : set.getMSnIdentificationNames()){
                String oneName = null;
                if (msnNames instanceof Vector){
                  oneName = ((Vector<String>)msnNames).get(0);
                }else{
                  oneName = ((String)msnNames);
                }
                if (!oneName.replaceAll("/","_").equalsIgnoreCase(ldaName))continue;
                if (msnNames instanceof Vector){
                  oneName = "";
                  for (String name: ((Vector<String>)msnNames)){
                    oneName+=name+";";
                  }
                  oneName = oneName.substring(0,oneName.length()-1);
                }else{
                  oneName = ((String)msnNames);
                }
                if (nameString==null){
                  nameDidNotChange = true;
                  nameString = oneName;
                }else if (nameDidNotChange){
                  if (!nameString.equalsIgnoreCase(oneName)){
                    nameDidNotChange = false;
                    nameString = ldaName;
                  }
                }
              }
              
              boolean isAtCorrectRt = isHitAtCorrectRt(Double.parseDouble(set.getRt()), info, ref);
              if (isAtCorrectRt){
                ms2AnalyteDetected = true;
                if (ldaMS2Evidence==-3 || ldaMS2Evidence==-2) ldaMS2Evidence=1;
                else if (ldaMS2Evidence==-1 || ldaMS2Evidence==0) ldaMS2Evidence=2;
              } else {
                if (ldaMS2Evidence==-2||ldaMS2Evidence==-1) ldaMS2Evidence=-3;
                else if (ldaMS2Evidence==0) ldaMS2Evidence=-2;
                else if (ldaMS2Evidence==2) ldaMS2Evidence=1;          
              }
              rtMS2LDA += set.getRt()+" ";
            }
            if (ldaMS2Evidence>0 || ldaMS2Evidence<-1) ldaMS2Name = nameString; 
          }
          if (rtMS2LDA.length()>0) rtMS2LDA.substring(0,rtMS2LDA.length()-1);
          
          // here starts the LB check of the MS2 evidence
          Hashtable<String,String> lbNames = new  Hashtable<String,String>();
          int lbMS2Evidence = -1;
          double lbMs2Prob = 0d;
          if (ident!=null){
            for (String lbName : ident.getIdentifications(info.getSns())){     
              if (!StaticUtils.isAPermutedVersion(keyWOSlash, lbName.replaceAll("/", "_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) continue;
              lbNames.put(lbName, lbName);
              lbOKAccordingToReference.put(lbName, lbName);
              double prob = ident.getHighestProbability(lbName, info.getSns());
              if (prob>lbMs2Prob) lbMs2Prob=prob;
              for (String rt : ident.getRetentionTimes(lbName, info.getSns())){
                boolean isAtCorrectRt = isHitAtCorrectRt(Double.parseDouble(rt), info, ref);
                if (isAtCorrectRt){
                  ms2AnalyteDetected = true;
                  if (lbMS2Evidence==-3 || lbMS2Evidence==-2) lbMS2Evidence=1;
                  else if (lbMS2Evidence==-1 || lbMS2Evidence==0) lbMS2Evidence=2;
                } else {
                  if (lbMS2Evidence==-2||lbMS2Evidence==-1) lbMS2Evidence=-3;
                  else if (lbMS2Evidence==0) lbMS2Evidence=-2;
                  else if (lbMS2Evidence==2) lbMS2Evidence=1;          
                }
                rtMS2LB += rt+" ";              
              }
            }
          }
          String lbName = "";
          if (lbMS2Evidence!=-1 && lbNames.size()>0){
            for (String name : lbNames.keySet()) lbName+=name+" ";
          }
          if (lbName.length()>1)lbName=lbName.substring(0,lbName.length()-1);
          else lbName = "not reported";
          if (rtMS2LB.length()>1) rtMS2LB = rtMS2LB.substring(0,rtMS2LB.length()-1);
          if (ldaMS2Evidence!=-1 || lbMS2Evidence!=-1){
            ms2Compare.add(new LdaLBLASTCompareVO(key,ldaMS2Name,lbName,ldaMS2Evidence,lbMS2Evidence,rtMS2Ideal,rtMS2LDA,rtMS2LB,false,
                lbMs2Prob, (float)quant.getAnalyteMass(), false));
          }
        }
      }
    }else{
      correctMS1Name += "_FP";
      for (LipidParameterSet set : ldaMs2SpectraIdentified){
        ldaMS1Name = ms1Analyte;
        ldaMS1Evidence = LdaLBLASTCompareVO.FALSE_POSITIVE;
      }
      LinkedHashMap<String,String> rts = new LinkedHashMap<String,String>();
      if (ident!=null && ident.getDetections().size()>0){
        lbMS1Name = ms1Analyte;
        lbMS1Evidence = LdaLBLASTCompareVO.FALSE_POSITIVE;
        for (LipidBLASTDetectionVO detect : ident.getDetections()){
          if (detect.getProbability()>lbProb) lbProb = detect.getProbability();
          rts.put(detect.getRetentionTime(), detect.getRetentionTime());
        }
      }
      if (rts.size()>0){
        for(String rt : rts.keySet()) rtLB += rt+" ";
        rtLB = rtLB.substring(0,rtLB.length()-1);
      }
      isFP = true;
    }
    //add now the ms2 FPs of the LDA including the ones that LipidBlast shows too
    if (info.getSns()>1 && (ldaMS1Evidence<-1 || ldaMS1Evidence>0)){
      for (String ldaKey : ldaMs2Evidence.keySet()){
        if (ldaOKAccordingToReference.containsKey(ldaKey)) continue;
        Vector<LipidomicsMSnSet> ldaHits = ldaMs2Evidence.get(ldaKey);
        String nameString = null;
        if (ldaHits.size()==0) continue;
        boolean nameDidNotChange = false;
        String rtMS2LDA = "";
        for (LipidomicsMSnSet set : ldaHits){
          for (Object msnNames : set.getMSnIdentificationNames()){
            String oneName = null;
            if (msnNames instanceof Vector){
              oneName = ((Vector<String>)msnNames).get(0);
            }else{
              oneName = ((String)msnNames);
            }
            if (!oneName.replaceAll("/","_").equalsIgnoreCase(ldaKey))continue;
            if (msnNames instanceof Vector){
              oneName = "";
              for (String name: ((Vector<String>)msnNames)){
                oneName+=name+";";
              }
              oneName = oneName.substring(0,oneName.length()-1);
            }else{
              oneName = ((String)msnNames);
            }
            if (nameString==null){
              nameDidNotChange = true;
              nameString = oneName;
            }else if (nameDidNotChange){
              if (!nameString.equalsIgnoreCase(oneName)){
                nameDidNotChange = false;
                nameString = ldaKey;
              }
            }
            rtMS2LDA += set.getRt()+" ";
          }
        }
        int lbMS2Evidence = 0;
        double lbMs2Prob = 0d;
        String rtMS2LB = "";
        Hashtable<String,String> lbNames = new  Hashtable<String,String>();
        //check if LipidBlast finds the FP too
        if (ident!=null && (lbMS1Evidence<-1  || lbMS1Evidence>0)){
          for (String lbNm : ident.getIdentifications(info.getSns())){
            if (!StaticUtils.isAPermutedVersion(ldaKey, lbNm.replaceAll("/", "_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) continue;
            lbMS2Evidence = LdaLBLASTCompareVO.FALSE_POSITIVE;
            lbNames.put(lbNm, lbNm);
            lbOKAccordingToReference.put(lbNm, lbNm);
            double prob = ident.getHighestProbability(lbNm, info.getSns());
            if (prob>lbMs2Prob) lbMs2Prob=prob;
            for (String rt : ident.getRetentionTimes(lbNm, info.getSns())){
              rtMS2LB += rt+" ";              
            }
          }         
        }
        String lbName = "";
        if (lbMS2Evidence!=-1 && lbNames.size()>0){
          for (String name : lbNames.keySet()) lbName+=name+" ";
        }
        if (lbName.length()>1)lbName=lbName.substring(0,lbName.length()-1);
        if (rtMS2LB.length()>1) rtMS2LB = rtMS2LB.substring(0,rtMS2LB.length()-1);
        ms2Compare.add(new LdaLBLASTCompareVO(nameString+"_FP",nameString,lbName,-2,lbMS2Evidence,"",rtMS2LDA,rtMS2LB,true,
            lbMs2Prob, (float)quant.getAnalyteMass(), false));
      }
    }
    //add now the ms2 FPs that were found by LipidBlast only
    if (info.getSns()>1 && ident!=null && (lbMS1Evidence<-1  || lbMS1Evidence>0)){
      int lbMS2Evidence = 0;
      Hashtable<String,String> lbNames = new  Hashtable<String,String>();
      Hashtable<String, Hashtable<String,String>> lbNameVariations = new Hashtable<String, Hashtable<String,String>>();
      Hashtable<String,String> rtMS2LBs = new Hashtable<String,String>();
      Hashtable<String,Double> probs = new Hashtable<String,Double>();
      for (String lbNm : ident.getIdentifications(info.getSns())){
        if (lbOKAccordingToReference.containsKey(lbNm)) continue;
        lbMS2Evidence = LdaLBLASTCompareVO.FALSE_POSITIVE;
        Hashtable<String,String> nameVariations = new Hashtable<String,String>();
        String rts = "";
        double prob = 0d;
        String key = lbNm;
        for (String other : lbNames.keySet()){
          if (StaticUtils.isAPermutedVersion(other.replaceAll("/", "_"),lbNm.replaceAll("/","_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
            nameVariations = lbNameVariations.get(other);
            rts = rtMS2LBs.get(other);
            key = other;
            break;
          }
        }
        nameVariations.put(lbNm, lbNm);
        double prb = ident.getHighestProbability(lbNm, info.getSns());
        if (prb>prob) prob = prb;
        for (String rt : ident.getRetentionTimes(lbNm, info.getSns())){
          rts += rt+" ";              
        }
        lbNames.put(key, lbNm);
        lbNameVariations.put(key, nameVariations);
        rtMS2LBs.put(key, rts);
        probs.put(key, prob);
      }
      if (lbMS2Evidence==LdaLBLASTCompareVO.FALSE_POSITIVE && lbNames.size()>0){
        for (String key : lbNames.keySet()){
          String lbName = "";
          Hashtable<String,String> nameVariations = lbNameVariations.get(key);
          String rtMS2LB = rtMS2LBs.get(key);
          double lbMs2Prob = probs.get(key);
          if (lbMS2Evidence!=-1 && lbNames.size()>0){
            for (String name : nameVariations.keySet()) lbName+=name+" ";
          }
          if (lbName.length()>1)lbName=lbName.substring(0,lbName.length()-1);
          if (rtMS2LB.length()>0) rtMS2LB = rtMS2LB.substring(0,rtMS2LB.length()-1);
          ms2Compare.add(new LdaLBLASTCompareVO(lbName+"_FP","not reported",lbName,0,lbMS2Evidence,"","",rtMS2LB,true,
            lbMs2Prob, (float)quant.getAnalyteMass(), false));
        }
      }
    }
    LdaLBLASTCompareVO compare = null;
    boolean createCompareVO = true;
    if (ldaMS1Evidence==0 && lbMS1Evidence==0 && !ms1Only) createCompareVO = false;
    if (ldaMS1Evidence<1 && lbMS1Evidence<1 && noMS2) correctMS1Name += "_noMS2";
    if (createCompareVO){
      compare = new LdaLBLASTCompareVO(correctMS1Name,ldaMS1Name,lbMS1Name,ldaMS1Evidence,lbMS1Evidence,rtIdeal,rtLDA,rtLB,isFP,
        lbProb, (float)quant.getAnalyteMass(), ms1Only);
      if (ms2Compare.size()>0) compare.setMs2Evidence(ms2Compare);
        
    }
    return compare;
    
/****    
    Vector<LdaLBLASTCompareVO> compareVOs = new Vector<LdaLBLASTCompareVO>();
    Vector<LdaLBLASTCompareVO> foundByLDAOnly = new Vector<LdaLBLASTCompareVO>();
    Vector<LdaLBLASTCompareVO> foundByLBOnly = new Vector<LdaLBLASTCompareVO>();
    //this was added for backward compatibility only, since now ms1 only hits are in the LDA results
    if (ldaMS1Only.size()==ldaAnalytes.size()&&ident==null) return compareVOs;
    
    Vector<LdaMinInfoVO> uniqueLDAInfos = new Vector<LdaMinInfoVO>();
    Hashtable<String,String> uniqueLDANames = new Hashtable<String,String>();
    boolean allMS1 = true;
    for (LipidParameterSet set : ldaAnalytes){
      if (!(set instanceof LipidomicsMSnSet)) continue;
      LipidomicsMSnSet ldaAnalyte = (LipidomicsMSnSet)set;
      for (Object msnNames : ldaAnalyte.getMSnIdentificationNames()){
        boolean ms1Only = ldaAnalyte.getStatus()==2 ? true : false;
        if (!ms1Only && allMS1) uniqueLDAInfos.clear();
        if (ms1Only && !allMS1) continue;
        String nameString = "";
        if (msnNames instanceof Vector){
          for (String name : (Vector<String>)msnNames){
            nameString += name+";";
          }
          nameString = nameString.substring(0,nameString.length()-1);
        }else{
          nameString = (String)msnNames;
        }
        if (uniqueLDANames.containsKey(nameString)){
          int posPartner = -1;
          for (int i=0;i!=uniqueLDAInfos.size();i++){
            if (uniqueLDAInfos.get(i).getIdentification().equalsIgnoreCase(nameString)){
              posPartner = i;
              break;
            }
          }
          if (ldaAnalyte.Area>uniqueLDAInfos.get(posPartner).getArea()){
            uniqueLDAInfos.remove(posPartner);
            uniqueLDAInfos.add(new LdaMinInfoVO(nameString,ldaAnalyte.getRt(),ldaAnalyte.Area,ms1Only));
            uniqueLDANames.put(nameString, nameString);
          }
        }else{
          uniqueLDAInfos.add(new LdaMinInfoVO(nameString,ldaAnalyte.getRt(),ldaAnalyte.Area,ms1Only));
          uniqueLDANames.put(nameString, nameString);
        }
        if (!ms1Only) allMS1 = false;
      } 
    }
    if (uniqueLDAInfos.size()>0){
      for (LdaMinInfoVO ldaVO : uniqueLDAInfos){
        String nameString = ldaVO.getIdentification();
        boolean foundInLBlast = false;
        String blastName = null;
        if (ident!=null){
          for (String found : ident.getIdentifications(info.getSns())){
            String toCompare = new String(nameString);
            if (toCompare.indexOf(";")!=-1) toCompare = toCompare.substring(0,toCompare.indexOf(";"));
            if (ldaVO.isMs1Only() || checkForFACorrectness(toCompare,found)[0]){
              foundInLBlast = true;
              blastName = found;
              break;
            }
          }
        }
        if (foundInLBlast) compareVOs.add(new LdaLBLASTCompareVO(nameString,blastName,ldaVO.getRt(),ldaVO.isMs1Only()));
        else foundByLDAOnly.add(new LdaLBLASTCompareVO(nameString,null,ldaVO.getRt(),ldaVO.isMs1Only()));
      }
      if (ident!=null){
      for (String found : ident.getIdentifications(info.getSns())){     
        boolean foundInLDA = false;
        for (LdaMinInfoVO ldaVO : uniqueLDAInfos){
          String nameString = ldaVO.getIdentification();
          String toCompare = new String(nameString);
          if (toCompare.indexOf(";")!=-1) toCompare = toCompare.substring(0,toCompare.indexOf(";"));
          if (ldaVO.isMs1Only() || checkForFACorrectness(toCompare,found)[0]){
            foundInLDA = true;
            break;
          } 
        }
        if (!foundInLDA){
          boolean alreadyAdded = false;
          for (LdaLBLASTCompareVO compVO : foundByLBOnly){
            if (checkForFACorrectness(compVO.getLbName(),found)[0]){
              alreadyAdded = false;
              break;
            }
          }
          if (!alreadyAdded) foundByLBOnly.add(new LdaLBLASTCompareVO(null,found,null,false));
        }
      }
      }
    } else {
      for (String found : ident.getIdentifications(info.getSns())){     
        boolean alreadyAdded = false;
        for (LdaLBLASTCompareVO compVO : foundByLBOnly){
          if (checkForFACorrectness(compVO.getLbName(),found)[0]){
            alreadyAdded = true;
            break;
          }
        }
        if (!alreadyAdded) foundByLBOnly.add(new LdaLBLASTCompareVO(null,found,null,false)); 
      }
    }
    setFullLBNamesInCompareVOs(compareVOs,ident,info.getSns());
    setFullLBNamesInCompareVOs(foundByLBOnly,ident,info.getSns());
    compareVOs.addAll(foundByLDAOnly);
    compareVOs.addAll(foundByLBOnly);
    return compareVOs;*/
  }

  private boolean isHitAtCorrectRt(double rt, LipidClassInfoVO info, LinkedHashMap<String,ReferenceInfoVO> species){
    boolean isAtCorrectRt = false;
    for (String ms2Name : species.keySet()){
      ReferenceInfoVO ref = species.get(ms2Name);
      if (isHitAtCorrectRt(rt,info,ref)){
        isAtCorrectRt = true;
        break;
      }
    }
    return isAtCorrectRt;
  }
  
  private boolean isHitAtCorrectRt(double rt, LipidClassInfoVO info, ReferenceInfoVO ref){
    if (!info.checkRt()) return true;
    boolean isAtCorrectRt = false;
    for (int i=0; i!=ref.getCorrectRts().length; i++) {
      double startRt = ref.getCorrectRts()[i]-info.getRtTolerance();
      double stopRt = ref.getCorrectRts()[i]+info.getRtTolerance();
      if (startRt<=rt && rt<=stopRt) isAtCorrectRt = true;
    }
    return isAtCorrectRt;    
  }

  
//  private Vector<LdaLBLASTCompareVO> getLdaLBlastSpeciesComparision(Vector<LipidParameterSet> ldaAnalytes, LipidBLASTIdentificationVO ident, int amountSns){
//    Vector<LdaLBLASTCompareVO> compareVOs = new Vector<LdaLBLASTCompareVO>();
//    Vector<LdaLBLASTCompareVO> foundByLDAOnly = new Vector<LdaLBLASTCompareVO>();
//    Vector<LdaLBLASTCompareVO> foundByLBOnly = new Vector<LdaLBLASTCompareVO>();
//    
//    Vector<LdaMinInfoVO> uniqueLDAInfos = new Vector<LdaMinInfoVO>();
//    Hashtable<String,String> uniqueLDANames = new Hashtable<String,String>();
//    boolean allMS1 = true;
//    for (LipidParameterSet set : ldaAnalytes){
//      LipidomicsMSnSet ldaAnalyte = (LipidomicsMSnSet)set;
//      for (Object msnNames : ldaAnalyte.getMSnIdentificationNames()){
//        boolean ms1Only = ldaAnalyte.getStatus()==2 ? true : false;
//        if (!ms1Only && allMS1) uniqueLDAInfos.clear();
//        if (ms1Only && !allMS1) continue;
//        String nameString = "";
//        if (msnNames instanceof Vector){
//          for (String name : (Vector<String>)msnNames){
//            nameString += name+";";
//          }
//          nameString = nameString.substring(0,nameString.length()-1);
//        }else{
//          nameString = (String)msnNames;
//        }
//        if (uniqueLDANames.containsKey(nameString)){
//          int posPartner = -1;
//          for (int i=0;i!=uniqueLDAInfos.size();i++){
//            if (uniqueLDAInfos.get(i).getIdentification().equalsIgnoreCase(nameString)){
//              posPartner = i;
//              break;
//            }
//          }
//          if (ldaAnalyte.Area>uniqueLDAInfos.get(posPartner).getArea()){
//            uniqueLDAInfos.remove(posPartner);
//            uniqueLDAInfos.add(new LdaMinInfoVO(nameString,ldaAnalyte.getRt(),ldaAnalyte.Area,ms1Only));
//            uniqueLDANames.put(nameString, nameString);
//          }
//        }else{
//          uniqueLDAInfos.add(new LdaMinInfoVO(nameString,ldaAnalyte.getRt(),ldaAnalyte.Area,ms1Only));
//          uniqueLDANames.put(nameString, nameString);
//        }
//        if (!ms1Only) allMS1 = false;
//      } 
//    }
//    if (uniqueLDAInfos.size()>0){
//      for (LdaMinInfoVO ldaVO : uniqueLDAInfos){
//        String nameString = ldaVO.getIdentification();
//        boolean foundInLBlast = false;
//        String blastName = null;
//        if (ident!=null){
//          for (String found : ident.getIdentifications(amountSns)){
//            String toCompare = new String(nameString);
//            if (toCompare.indexOf(";")!=-1) toCompare = toCompare.substring(0,toCompare.indexOf(";"));
//            if (ldaVO.isMs1Only() || checkForFACorrectness(toCompare,found)[0]){
//              foundInLBlast = true;
//              blastName = found;
//              break;
//            }
//          }
//        }
//        if (foundInLBlast) compareVOs.add(new LdaLBLASTCompareVO(nameString,blastName,ldaVO.getRt(),ldaVO.isMs1Only()));
//        else foundByLDAOnly.add(new LdaLBLASTCompareVO(nameString,null,ldaVO.getRt(),ldaVO.isMs1Only()));
//      }
//      if (ident!=null){
//      for (String found : ident.getIdentifications(amountSns)){     
//        boolean foundInLDA = false;
//        for (LdaMinInfoVO ldaVO : uniqueLDAInfos){
//          String nameString = ldaVO.getIdentification();
//          String toCompare = new String(nameString);
//          if (toCompare.indexOf(";")!=-1) toCompare = toCompare.substring(0,toCompare.indexOf(";"));
//          if (ldaVO.isMs1Only() || checkForFACorrectness(toCompare,found)[0]){
//            foundInLDA = true;
//            break;
//          } 
//        }
//        if (!foundInLDA){
//          boolean alreadyAdded = false;
//          for (LdaLBLASTCompareVO compVO : foundByLBOnly){
//            if (checkForFACorrectness(compVO.getLbName(),found)[0]){
//              alreadyAdded = false;
//              break;
//            }
//          }
//          if (!alreadyAdded) foundByLBOnly.add(new LdaLBLASTCompareVO(null,found,null,false));
//        }
//      }
//      }
//    } else {
//      for (String found : ident.getIdentifications(amountSns)){     
//        boolean alreadyAdded = false;
//        for (LdaLBLASTCompareVO compVO : foundByLBOnly){
//          if (checkForFACorrectness(compVO.getLbName(),found)[0]){
//            alreadyAdded = true;
//            break;
//          }
//        }
//        if (!alreadyAdded) foundByLBOnly.add(new LdaLBLASTCompareVO(null,found,null,false)); 
//      }
//    }
//    setFullLBNamesInCompareVOs(compareVOs,ident,amountSns);
//    setFullLBNamesInCompareVOs(foundByLBOnly,ident,amountSns);
//    compareVOs.addAll(foundByLDAOnly);
//    compareVOs.addAll(foundByLBOnly);
//    return compareVOs;
//  }
/****  
  private void setFullLBNamesInCompareVOs(Vector<LdaLBLASTCompareVO> comps, LipidBLASTIdentificationVO ident, int amountSns){
    for (LdaLBLASTCompareVO comp : comps){
      String finalName = "";
      String finalRt = "";
      Hashtable<String,String> alreadyAdded = new Hashtable<String,String>();
      Hashtable<String,String> rtAdded = new Hashtable<String,String>();
      double highestProb = 0d;
      for (String name : ident.getIdentifications(amountSns)){
        if ((comp.isLdaMs1Only_() || checkForFACorrectness(comp.getLbName(),name)[0]) && !alreadyAdded.containsKey(name)){          
          finalName += name+" ";
          alreadyAdded.put(name, name);
          if (ident.getHighestProbability(name, amountSns)>highestProb) highestProb = ident.getHighestProbability(name, amountSns);
          for (String rt : ident.getRetentionTimes(name, amountSns)){
            if (!rtAdded.containsKey(rt)){
              finalRt += rt+" ";
              rtAdded.put(rt, rt);
            }
          }
        }
      }
      if (finalName.length()>0){
        finalName = finalName.substring(0,finalName.length()-1);
        finalRt = finalRt.substring(0,finalRt.length()-1);
        comp.setLbName(finalName);
        comp.setLbProb(highestProb);
        comp.setLbRts(finalRt);
        if (alreadyAdded.size()>1) comp.setLbIsAmbiguous(true);
      }
    }
  }*/
  
  //TODO: is not corrected yet
  private void checkForUniqueFACombinations(Vector<LdaLBLASTCompareVO> compares){
    Hashtable<String,Boolean> fasUsed = new Hashtable<String,Boolean>(); 
    
    Vector<LdaLBLASTCompareVO> left = new Vector<LdaLBLASTCompareVO>();
    for (LdaLBLASTCompareVO compare : compares){
      if (!compare.isLdaFound()) continue;
      String[] fas = compare.getLdaName().replaceAll("/", "_").split("_");
      for (String fa:fas) fasUsed.put(fa, false);
      left.add(compare);
    }
    int count=0;
    while (left.size()>0 /*&& count<1*/){
      Vector<Integer> toRemove = new Vector<Integer>();
      int notUsed = 0;
      for (boolean value : fasUsed.values()) if (!value)notUsed++;
      //check if all of the fatty acids were used already
      for (int i=0;i!=left.size();i++){
        LdaLBLASTCompareVO compare = left.get(i);
        String[] fas = compare.getLdaName().replaceAll("/", "_").split("_");
        boolean allUsed = true;
        for (String fa:fas){
          if (!fasUsed.get(fa)) allUsed=false;
        }
        if (allUsed){
          compare.setFasInOtherCombination(true);
          toRemove.add(i);
        }
      }
      for (int i=(toRemove.size()-1);i!=-1;i--) left.removeElementAt(toRemove.get(i));
      
      //build statistics on FA usage
      Hashtable<String,Integer> amountSameFAUsedInCombination = new Hashtable<String,Integer>();
      for (int i=0;i!=left.size();i++){
        LdaLBLASTCompareVO compare = left.get(i);
        String[] fas = compare.getLdaName().replaceAll("/", "_").split("_");
        Hashtable<String,String> usedFA = new Hashtable<String,String>();
        for (String fa:fas){
          if (usedFA.containsKey(fa)) continue;
          if (fasUsed.get(fa)) continue;
          int amount = 0;
          if (amountSameFAUsedInCombination.containsKey(fa)) amount = amountSameFAUsedInCombination.get(fa);
          amount++;
          amountSameFAUsedInCombination.put(fa,amount);
          usedFA.put(fa, fa);
        }
      }
      boolean areThereSingleFAs = false;
      toRemove = new Vector<Integer>();
      Hashtable<Integer,Integer> used = new Hashtable<Integer,Integer>();
      for (String fa : amountSameFAUsedInCombination.keySet()){
        if (amountSameFAUsedInCombination.get(fa)==1){
          areThereSingleFAs = true;
          for (int i=0;i!=left.size();i++){
            LdaLBLASTCompareVO compare = left.get(i);
            String[] fas = compare.getLdaName().replaceAll("/", "_").split("_");
            boolean containsFA = false;
            for (String otherFA:fas){
              if (fa.equalsIgnoreCase(otherFA)){
                containsFA = true;
                break;
              }
            }
            if (!containsFA) continue;
            if (used.containsKey(i))continue;
            toRemove.add(i);
            used.put(i, i);
            for (String otherFA:fas){
              fasUsed.put(otherFA, true);
            }
          }
        }
      }
      Collections.sort(toRemove);
      for (int i=(toRemove.size()-1);i!=-1;i--)left.removeElementAt(toRemove.get(i));
      
      // if there are no single FAS, we have to decide for one FA hit
      if (!areThereSingleFAs && left.size()>0 && amountSameFAUsedInCombination.size()>0){
        int lowestFANumber = Integer.MAX_VALUE;
        Hashtable<String,String> faCandidates = new Hashtable<String,String>();
        for (String fa : amountSameFAUsedInCombination.keySet()){
          int foundIn = amountSameFAUsedInCombination.get(fa);
          if (foundIn<lowestFANumber){
            lowestFANumber = foundIn;
            faCandidates = new Hashtable<String,String>();
            faCandidates.put(fa,fa);
          } else if (foundIn<lowestFANumber) faCandidates.put(fa,fa);
        }
        //check which combinations contain such a candidate, and which has the
        //highest number of diverse FAs
        int numberOfDiverseFAs = 0;
        int oneToRemove = -1;
        for (int i=0;i!=left.size();i++){
          LdaLBLASTCompareVO compare = left.get(i);
          String[] fas = compare.getLdaName().replaceAll("/", "_").split("_");
          boolean containsFA = false;
          for (String otherFA:fas){
            if (faCandidates.containsKey(otherFA)){
              containsFA = true;
              break;
            }
          }
          if (containsFA){
            Hashtable<String,String> usedFA = new Hashtable<String,String>();
            int unusedFAs = 0;
            for (String otherFA:fas){
              if (usedFA.containsKey(otherFA)) continue;
              usedFA.put(otherFA, otherFA);
              if (!fasUsed.get(otherFA)) unusedFAs++;
            }
            if (unusedFAs>numberOfDiverseFAs){
              numberOfDiverseFAs = unusedFAs;
              oneToRemove = i;
            }
          }
        }
        String[] fas = left.get(oneToRemove).getLdaName().replaceAll("/", "_").split("_");
        for (String otherFA:fas){
          fasUsed.put(otherFA, true);
        }
        left.removeElementAt(oneToRemove);
      }
      count++;
    }
  }
  
  private void testAfterStructuralIdentification(){
    try{
      String[] chromPaths = StringUtils.getChromFilePaths("D:\\ABSciex\\20150501\\D_2.chrom");   
      LipidomicsAnalyzer analyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);

      Hashtable<String,Vector<LipidParameterSet>> idents = LDAResultReader.readResultFile("D:\\ABSciex\\20150501\\D_2_TG_DG_PC_PE_SM_sent.xlsx", new Hashtable<String,Boolean>()).getIdentifications();
      long time = System.currentTimeMillis();
      for (String className : idents.keySet()){
        for (LipidParameterSet set : idents.get(className)){
//          System.out.println(set.getNameStringWithoutRt());
          MSnAnalyzer msnAnalyzer = new MSnAnalyzer(className,set.getModificationName(),set,analyzer,null,true,false);  
          if (msnAnalyzer.checkStatus()!=LipidomicsMSnSet.DISCARD_HIT) msnAnalyzer.getResult();

        }
      }
      long timePerThread = (System.currentTimeMillis()-time)/7l;
      System.out.println("Required time: "+(timePerThread/(60*1000))+" minutes "+timePerThread%(60*1000)/1000+" seconds");

      
    } catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void checkForMassDeviation(){
    String basePath = "D:\\Experiment1\\QTOF\\positive\\";
    File[] files = (new File(basePath)).listFiles();
    files = new File[1];
    //files[0] = new File("D:\\Experiment1\\SanDiego\\positive\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Standard 1 5uM pos_Ex1_pos_new.xlsx");
    files[0] = new File("D:\\Experiment1\\SanDiego\\negative\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Standard 1 5uM neg_Ex1_neg.xlsx");

    try {
      ////Vector allInfo = QuantificationThread.parseQuantExcelFile("D:\\Experiment1\\massLists\\quant\\Ex1_pos.xlsx",  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f);
      Vector allInfo = null;
      Set<String> classes = ((LinkedHashMap<String,Integer>)allInfo.get(0)).keySet();
      Hashtable<String,Vector<String>> analyteSequence = (Hashtable<String,Vector<String>>)allInfo.get(1);
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>) allInfo.get(3);

      Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>> results = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>>();
      for (int i=0; i!=files.length; i++){
        File file = files[i];
        if (!file.isFile()) continue;
        if (!file.getName().endsWith(".xlsx")) continue;
        Hashtable<String,Vector<LipidParameterSet>> idents = LDAResultReader.readResultFile(file.getAbsolutePath(), new Hashtable<String,Boolean>()).getIdentifications();
        for (String className : idents.keySet()){
          Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> analytes = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>>();
          if (results.containsKey(className)) analytes = results.get(className);
          Vector<LipidParameterSet> sets = idents.get(className);
          for (LipidParameterSet set : sets){
            if (isRtInInterestingRegion(className, set)){
              String id = set.getName()+":"+set.getDoubleBonds().toString();
              Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> mods = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
              if (analytes.containsKey(id)) mods = analytes.get(id);
              Hashtable<String,Hashtable<String,LipidParameterSet>> fileResults = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
              if (mods.containsKey(set.getModificationName())) fileResults = mods.get(set.getModificationName());
              Hashtable<String,LipidParameterSet> rts = new Hashtable<String,LipidParameterSet>();
              if (fileResults.containsKey(file.getName())) rts = fileResults.get(file.getName());
              rts.put(set.getRt(), set);
              fileResults.put(file.getName(),rts);
              mods.put(set.getModificationName(),fileResults);
              analytes.put(id, mods);
            }
          }
          results.put(className, analytes);
        }
      }
      for (String className: classes){
        Vector<String> analytesInSequence = analyteSequence.get(className);
        Hashtable<String,Hashtable<String,QuantVO>> analyteQuants = quantObjects.get(className);
        Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>> analyteRes = results.get(className);
        System.out.println("!!!!!! "+className+" !!!!!!");
        for (String anal : analytesInSequence){
          Hashtable<String,QuantVO> modQuants = analyteQuants.get(anal);
          if (analyteRes==null || !analyteRes.containsKey(anal)) continue;
          Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> modRes = analyteRes.get(anal);
          for (String mod : modQuants.keySet()){
            QuantVO quant = modQuants.get(mod);
            if (!modRes.containsKey(mod)) continue;
            Hashtable<String,Hashtable<String,LipidParameterSet>> fileResults = modRes.get(mod);
            for (String fileName : fileResults.keySet()){
              Hashtable<String,LipidParameterSet> rts = fileResults.get(fileName);
              for (String rt : rts.keySet()){
                LipidParameterSet set = rts.get(rt);
                Vector<CgProbe> probes = set.getIsotopicProbes().get(0);
                for (CgProbe probe : probes){
                  float difference = (probe.Mz-(float)quant.getAnalyteMass())*1000f;
                  if (difference<0f) difference = difference*-1f;
                  if (difference>150f)
                    System.out.println(fileName+": "+anal+"_"+mod+"_"+rt+": "+(probe.Mz-(float)quant.getAnalyteMass())*1000f);
                }
              }
            }
          }
        }
      }
    }
    catch (Exception e) {
      e.printStackTrace();
    }
    
  }
  
  private boolean isRtInInterestingRegion(String className, LipidParameterSet set){
    float rt = Float.parseFloat(set.getRt());
    String id = set.getName()+":"+set.getDoubleBonds().toString();
    if (className.equalsIgnoreCase("TG")){
      if (id.equalsIgnoreCase("44:1") && 27.7f<rt && rt<28.7f) return true;
      else if (id.equalsIgnoreCase("48:0") && 39f<rt && rt<40f) return true;
      else if (id.equalsIgnoreCase("48:1") && 33.4f<rt && rt<34.4f) return true;
      else if (id.equalsIgnoreCase("50:0") && 45.5f<rt && rt<46.5f) return true;
      else if (id.equalsIgnoreCase("50:1") && 38.1f<rt && rt<39.1f) return true;
      else if (id.equalsIgnoreCase("51:1") && 40.7f<rt && rt<41.7f) return true;
      else if (id.equalsIgnoreCase("58:7") && 29f<rt && rt<30f) return true;
      else if (id.equalsIgnoreCase("58:10") && 25.2f<rt && rt<26.2f) return true;
      else if (id.equalsIgnoreCase("60:1") && 39f<rt && rt<40f) return true;
      else if (id.equalsIgnoreCase("62:16") && 22.5f<rt && rt<23.5f) return true;
    } else if (className.equalsIgnoreCase("PI")){
      if (id.equalsIgnoreCase("25:0") && 18.9f<rt && rt<19.9f) return true;
      else if (id.equalsIgnoreCase("31:1") && 20.5f<rt && rt<21.5f) return true;
      else if (id.equalsIgnoreCase("37:4") && 21f<rt && rt<22f) return true;
      else if (id.equalsIgnoreCase("43:6") && 22f<rt && rt<23f) return true;
    } else if (className.equalsIgnoreCase("P-PC")){
      if (id.equalsIgnoreCase("36:1") && 24.1f<rt && rt<25.1f) return true;
      else if (id.equalsIgnoreCase("38:4") && 22.2f<rt && rt<23.2f) return true;
      else if (id.equalsIgnoreCase("40:6") && 21.9f<rt && rt<22.9f) return true;    
    } else if (className.equalsIgnoreCase("P-PE")){
      if (id.equalsIgnoreCase("36:1") && 24.2f<rt && rt<25.2f) return true;
      else if (id.equalsIgnoreCase("38:4") && 22.5f<rt && rt<23.5f) return true;
      else if (id.equalsIgnoreCase("40:6") && 22.0f<rt && rt<23.0f) return true;    
    } else if (className.equalsIgnoreCase("LPC")){
      if (id.equalsIgnoreCase("13:0") && 13.0f<rt && rt<14.0f) return true;
      else if (id.equalsIgnoreCase("14:0") && 14.4f<rt && rt<15.4f) return true;
      else if (id.equalsIgnoreCase("16:0") && 16.3f<rt && rt<17.3f) return true;
      else if (id.equalsIgnoreCase("18:0") && 17.6f<rt && rt<18.6f) return true;
      else if (id.equalsIgnoreCase("18:1") && 16.5f<rt && rt<17.5f) return true;
    } else if (className.equalsIgnoreCase("LPE")){
      if (id.equalsIgnoreCase("14:0") && 14.5f<rt && rt<15.5f) return true;
      else if (id.equalsIgnoreCase("18:0") && 17.8f<rt && rt<18.8f) return true;
      else if (id.equalsIgnoreCase("18:1") && 16.7f<rt && rt<17.7f) return true;
    } else if (className.equalsIgnoreCase("LPE")){
      if (id.equalsIgnoreCase("17:1") && 16.9f<rt && rt<17.9f) return true;
      else if (id.equalsIgnoreCase("18:0") && 17.1f<rt && rt<18.1f) return true;
      else if (id.equalsIgnoreCase("18:1") && 16.1f<rt && rt<17.1f) return true;
    } else if (className.equalsIgnoreCase("PC")){
      if (id.equalsIgnoreCase("32:0") && 22.3f<rt && rt<23.3f) return true;
      else if (id.equalsIgnoreCase("32:1") && 21.5f<rt && rt<22.5f) return true;
      else if (id.equalsIgnoreCase("34:0") && 23.8f<rt && rt<24.8f) return true;
      else if (id.equalsIgnoreCase("34:1") && 22.4f<rt && rt<23.3f) return true;
      else if (id.equalsIgnoreCase("36:0") && 25.2f<rt && rt<26.2f) return true;
      else if (id.equalsIgnoreCase("36:1") && 23.4f<rt && rt<24.4f) return true;
      else if (id.equalsIgnoreCase("36:2") && 22.3f<rt && rt<23.3f) return true;
      else if (id.equalsIgnoreCase("40:0") && 29.0f<rt && rt<30.9f) return true;
      else if (id.equalsIgnoreCase("48:2") && 34f<rt && rt<37f) return true;
    } else if (className.equalsIgnoreCase("PE")){
      if (id.equalsIgnoreCase("32:0") && 22.5f<rt && rt<23.5f) return true;
      else if (id.equalsIgnoreCase("34:0") && 23.7f<rt && rt<24.7f) return true;
      else if (id.equalsIgnoreCase("34:1") && 22.4f<rt && rt<23.4f) return true;
      else if (id.equalsIgnoreCase("36:0") && 25.1f<rt && rt<26.1f) return true;
      else if (id.equalsIgnoreCase("36:1") && 23.5f<rt && rt<24.5f) return true;
      else if (id.equalsIgnoreCase("36:2") && 22.5f<rt && rt<23.5f) return true;
      else if (id.equalsIgnoreCase("36:4") && 21.3f<rt && rt<22.3f) return true;
    } else if (className.equalsIgnoreCase("PG")){
      if (id.equalsIgnoreCase("32:0") && 21.4f<rt && rt<22.4f) return true;
      else if (id.equalsIgnoreCase("34:0") && 22.2f<rt && rt<23.2f) return true;
      else if (id.equalsIgnoreCase("34:1") && 21.3f<rt && rt<22.3f) return true;
      else if (id.equalsIgnoreCase("36:0") && 23.0f<rt && rt<24.0f) return true;
      else if (id.equalsIgnoreCase("36:1") && 22.1f<rt && rt<23.1f) return true;
      else if (id.equalsIgnoreCase("36:2") && 21.5f<rt && rt<22.5f) return true;
    } else if (className.equalsIgnoreCase("DG")){
//      if (id.equalsIgnoreCase("24:0") && 13f<rt && rt<15f) return true;
//      else if (id.equalsIgnoreCase("32:0") && 21.3f<rt && rt<22.3f) return true;
//      else if (id.equalsIgnoreCase("34:0") && 22.9f<rt && rt<23.9f) return true;
//      else if (id.equalsIgnoreCase("34:1") && 19.6f<rt && rt<20.6f) return true;
//      else if (id.equalsIgnoreCase("36:0") && 24.5f<rt && rt<25.5f) return true;
//      else if (id.equalsIgnoreCase("36:2") && 22f<rt && rt<22.5f) return true;
//      else if (id.equalsIgnoreCase("36:4") && 19.8f<rt && rt<20.8f) return true;
//      else if (id.equalsIgnoreCase("38:0") && 26f<rt && rt<27f) return true;
    } else if (className.equalsIgnoreCase("SM")){
      if (id.equalsIgnoreCase("18:0") && 22.9f<rt && rt<23.9f) return true;
      else if (id.equalsIgnoreCase("18:1") && 21.5f<rt && rt<22.5f) return true;
      else if (id.equalsIgnoreCase("24:0") && 28.5f<rt && rt<29.5f) return true;
      else if (id.equalsIgnoreCase("24:1") && 25.6f<rt && rt<26.6f) return true;
    } else if (className.equalsIgnoreCase("Cer")){
      if (id.equalsIgnoreCase("16:0") && 21.8f<rt && rt<22.8f) return true;
      else if (id.equalsIgnoreCase("17:0") && 22.2f<rt && rt<23.2f) return true;
      else if (id.equalsIgnoreCase("18:0") && 22.6f<rt && rt<23.6f) return true;
      else if (id.equalsIgnoreCase("18:1") && 21.8f<rt && rt<22.8f) return true;
      else if (id.equalsIgnoreCase("20:0") && 23.7f<rt && rt<24.8f) return true;
      
    }

    return false;
  }
    
  private void parseRule(){
    try{
      ElementConfigParser elPar = Settings.getElementParser();
      FragRuleParser parser = new FragRuleParser(elPar);
      parser.parseFile(new File("fragRules/PG_Na.frag.txt"));
    } catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void countMS2(){
    String filePath = "D:\\BiologicalExperiment\\Orbitrap_CID\\positive\\002_liver2-1_Orbitrap_CID_pos_TG_45_new.xlsx";
    Hashtable<String,String> ms2Idents= new Hashtable<String,String>();
    try {
      Hashtable<String,Boolean> showMods = new  Hashtable<String,Boolean>();
      QuantificationResult result = LDAResultReader.readResultFile(filePath, showMods);
      for (String key : result.getIdentifications().keySet()){
        if (!key.equalsIgnoreCase("TG")) continue;
        for (LipidParameterSet set : result.getIdentifications().get(key)){
          if (!(set instanceof LipidomicsMSnSet)) continue;
          LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
          Vector<Object> identNames = msn.getMSnIdentificationNames();
          for (Object nameObject : identNames){
            if (nameObject instanceof Vector){
              for (String name : (Vector<String>)nameObject){
                if (!ms2Idents.containsKey(name)) ms2Idents.put(name,name);
              }
            }else{
              String name = (String)nameObject;
              if (!ms2Idents.containsKey(name)) ms2Idents.put(name,name);
            }
          }
        }
      }
    }catch(Exception ex){
      ex.printStackTrace();
    }
    System.out.println(ms2Idents.size());
  }
  
  private void mergeTGRessults(){
    String basePath = "D:\\BiologicalExperiment\\Orbitrap_CID\\positive\\";
    Vector<String> files = new Vector<String>();
    files.add("002_liver2-1_Orbitrap_CID_pos_TG.xlsx");
    files.add("003_liver2-1_Orbitrap_CID_pos_TG.xlsx");
    files.add("004_liver2-1_Orbitrap_CID_pos_TG.xlsx");
    files.add("005_liver2-1_Orbitrap_CID_pos_TG.xlsx");
    files.add("006_liver2-1_Orbitrap_CID_pos_TG.xlsx");
    //read Result files and check how many different hits with different RTs are present (-> influences the amount of columns) + order the data + record the intensities of the msn hits
    Hashtable<String,Vector<LipidParameterSet>> allResults = new Hashtable<String,Vector<LipidParameterSet>>();
    Hashtable<String,Integer> amountDiffRts = new Hashtable<String,Integer>();
    // first key is the modification, second key is the ms1 name, third key is the file Name, fourth key is the retention time
    Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidomicsMSnSet>>>> orderedHits = new Hashtable<String,Hashtable<String,Hashtable<String,Hashtable<String,LipidomicsMSnSet>>>>();
    // first key is the modification, second key is the ms1 name, third key is ms2 name, value is the total area
    Hashtable<String,Hashtable<String,Hashtable<String,Double>>> ms2Areas = new Hashtable<String,Hashtable<String,Hashtable<String,Double>>>();
    for (String fileName : files){
      //first key is the modification name, the second key the species name; the third key the retention time
      int maxRts = 0;
      Hashtable<String,Hashtable<String,Hashtable<String,String>>> rtsPerModAndSpecies = new Hashtable<String,Hashtable<String,Hashtable<String,String>>>();
      try {
        Hashtable<String,Boolean> showMods = new  Hashtable<String,Boolean>();
        Vector<LipidParameterSet> result = LDAResultReader.readResultFile(basePath+fileName, showMods).getIdentifications().get("TG");
        allResults.put(fileName, result);
        for (LipidParameterSet set : result){
          if (!(set instanceof LipidomicsMSnSet)) continue;
          LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
          
          //check for the amount of different RTs 
          Hashtable<String,Hashtable<String,String>> speciesRTs = new Hashtable<String,Hashtable<String,String>>();
          if (rtsPerModAndSpecies.containsKey(msn.getModificationFormula())) speciesRTs = rtsPerModAndSpecies.get(msn.getModificationFormula());
          Hashtable<String,String> rts = new Hashtable<String,String>();
          if (speciesRTs.containsKey(set.getNameStringWithoutRt())) rts = speciesRTs.get(set.getNameStringWithoutRt());
          rts.put(set.getRt(), set.getRt());
          if (rts.size()>maxRts) maxRts = rts.size();
          speciesRTs.put(set.getNameStringWithoutRt(),rts);
          rtsPerModAndSpecies.put(msn.getModificationFormula(), speciesRTs);
          
          //put the hits in correct order
          Hashtable<String,Hashtable<String,Hashtable<String,LipidomicsMSnSet>>> oneMod = new Hashtable<String,Hashtable<String,Hashtable<String,LipidomicsMSnSet>>>();
          if (orderedHits.containsKey(msn.getModificationName())) oneMod = orderedHits.get(msn.getModificationName());
          Hashtable<String,Hashtable<String,LipidomicsMSnSet>> hitsOfSpecies = new Hashtable<String,Hashtable<String,LipidomicsMSnSet>>();
          if (oneMod.containsKey(msn.getNameStringWithoutRt())) hitsOfSpecies = oneMod.get(msn.getNameStringWithoutRt());
          Hashtable<String,LipidomicsMSnSet> hitsOfFile = new Hashtable<String,LipidomicsMSnSet>();
          if (hitsOfSpecies.containsKey(fileName)) hitsOfFile = hitsOfSpecies.get(fileName);
          hitsOfFile.put(msn.getRt(), msn);
          hitsOfSpecies.put(fileName, hitsOfFile);
          oneMod.put(msn.getNameStringWithoutRt(), hitsOfSpecies);
          orderedHits.put(msn.getModificationName(), oneMod);
          
          //get the Areas
          Hashtable<String,Hashtable<String,Double>> areasOneMod = new Hashtable<String,Hashtable<String,Double>>();
          if (ms2Areas.containsKey(msn.getModificationName())) areasOneMod = ms2Areas.get(msn.getModificationName());
          Hashtable<String,Double> areasMS1Species = new Hashtable<String,Double>();
          if (areasOneMod.containsKey(msn.getNameStringWithoutRt())) areasMS1Species = areasOneMod.get(msn.getNameStringWithoutRt());
          for (Object nameObject : msn.getMSnIdentificationNames()){
            String identificationString = "";
            double area = 0d;
            boolean noChain = false;
            if (nameObject instanceof Vector){
              area = msn.getRelativeIntensity(((Vector<String>)nameObject).get(0))*((double)msn.Area);
              for (String name : (Vector<String>)nameObject){
                identificationString+=name+";";
              }
              identificationString = identificationString.substring(0,identificationString.length()-1);
            }else{
              String name = (String) nameObject;
              identificationString = name;
              if (msn.getStatus()==LipidomicsMSnSet.HEAD_GROUP_DETECTED) noChain = true;
              else area = msn.getRelativeIntensity(name)*((double)msn.Area);
            }
            if (noChain) continue;
            if (areasMS1Species.containsKey(identificationString)) area+= areasMS1Species.get(identificationString);
            areasMS1Species.put(identificationString, area);
          }
          areasOneMod.put(msn.getNameStringWithoutRt(), areasMS1Species);
          ms2Areas.put(msn.getModificationName(), areasOneMod);
        }
      }catch(Exception ex){
        ex.printStackTrace();
      }
      amountDiffRts.put(fileName, maxRts);
    }
    int column = 3;
    Hashtable<String,Integer> startColumn = new Hashtable<String,Integer>();
    for (String file : files){
      startColumn.put(file, column);
      column += amountDiffRts.get(file)*2;
    }
    
    try{
      Vector<String> tgSpecies = null;
      ////tgSpecies = ((Hashtable<String,Vector<String>>)QuantificationThread.parseQuantExcelFile("D:\\BiologicalExperiment\\massLists\\positive\\TG.xlsx",  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f).get(1)).get("TG");
      Workbook workbook = null;
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(basePath+"Orbitrap_CID_TG_species_MSMS.xlsx"));
      Workbook wb = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(wb);
      CellStyle centerStyle = getCenterStyle(wb);
      Sheet sheet = wb.createSheet("TG_MS2");
      Vector<String> modifications = new Vector<String>();
      modifications.add("NH4");
      modifications.add("Na");
      int rowCount = 0;
      int count = 0;
      for (String mod : modifications){
        Row row = sheet.createRow(rowCount);
        rowCount++;
        count = 0;
        Cell cell = row.createCell(count);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("TG_"+mod+"_MS2");
        for (String fileName : files){
          count = startColumn.get(fileName);
          for (int i=0; i!=amountDiffRts.get(fileName); i++){
            count++;
            cell = row.createCell(count);
            String name = fileName.substring(0,fileName.indexOf("_Orbitrap_CID_pos_TG.xlsx"));
            cell.setCellValue(name+"_"+String.valueOf(i+1));
            count++;
          }
        }
        
        Hashtable<String,Hashtable<String,Hashtable<String,LipidomicsMSnSet>>> hits = orderedHits.get(mod);
        Hashtable<String,Hashtable<String,Double>> areas = ms2Areas.get(mod);
        for (String tgMS1 : tgSpecies){
          if (!hits.containsKey(tgMS1)) continue;
          Hashtable<String,Hashtable<String,LipidomicsMSnSet>> hitsOfMS1Species = hits.get(tgMS1);
          Hashtable<String,Double> ms2Hits = areas.get(tgMS1);
          Vector<String> ms2InDescendingOrder = new Vector<String>();
          for (String ms2 : ms2Hits.keySet()){
            int addPosition = 0;
            for (String other : ms2InDescendingOrder){
              if (ms2Hits.get(ms2)>ms2Hits.get(other)) break;
              addPosition++;
            }
            ms2InDescendingOrder.add(addPosition, ms2);
          }
          for (int i=0; i!=ms2InDescendingOrder.size(); i++){
            String ms2 = ms2InDescendingOrder.get(i);
            row = sheet.createRow(rowCount);
            rowCount++;
            count = 0;
            if (i==0){
              cell = row.createCell(count);
              cell.setCellValue(tgMS1);
            }
            count++;
            count++;
            cell = row.createCell(count);
            cell.setCellValue(ms2);
            for (String fileName : files){
              count = startColumn.get(fileName);
              if (hitsOfMS1Species.containsKey(fileName)){
                Hashtable<String,LipidomicsMSnSet> rtHits = hitsOfMS1Species.get(fileName);
                Vector<LipidomicsMSnSet> hitsInAscendingRt = new Vector<LipidomicsMSnSet>();
                for (String rt : rtHits.keySet()){
                  int add = 0;
                  for (LipidomicsMSnSet set : hitsInAscendingRt){
                    if (Float.parseFloat(rt)<Float.parseFloat(set.getRt()))break;
                    add++;
                  }
                  hitsInAscendingRt.add(add,rtHits.get(rt));
                }
                for (LipidomicsMSnSet set : hitsInAscendingRt){
                  Hashtable<String,Double> foundCombis = set.getChainCombinationRelativeAreas();
                  if (foundCombis.containsKey(ms2.replaceAll("/", "_"))){
                    cell = row.createCell(count);
                    cell.setCellValue(set.getRt());
                    count++;
                    cell = row.createCell(count);
                    cell.setCellStyle(centerStyle);
                    cell.setCellValue("found");
                    count++;                  
                  }else{
                    count++;
                    count++;
                  }
                }
              }else{
                count++;
                cell = row.createCell(count);
                cell.setCellValue("noMS2");
                cell.setCellStyle(centerStyle);
              }
            }
          }
        }
        
        rowCount++;
        rowCount++;
      }
      wb.write(out);
      out.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void setIgnoresForSpeciesLevelEvaluation(Vector<LdaLBLASTCompareVO> compares, boolean checkMS2){
    if (compares.size()<2) return;
    int ldaHighestEvidence = -4;
    int ldaHighestHit = 0;
    int lbHighestEvidence = -4;
    int lbHighestHit = 0;
    Hashtable<String,Integer> ldaMS2Ev = new Hashtable<String,Integer>();
    Hashtable<String,Integer> ldaMS2Hit = new Hashtable<String,Integer>();
    Hashtable<String,Integer> lbMS2Ev = new Hashtable<String,Integer>();
    Hashtable<String,Integer> lbMS2Hit = new Hashtable<String,Integer>();
    for (int i=0;i!=compares.size();i++){
      LdaLBLASTCompareVO comp = compares.get(i);
      if (isEvidenceStronger(ldaHighestEvidence,comp.getLdaIdentCode())){
        ldaHighestEvidence = comp.getLdaIdentCode();
        ldaHighestHit = i;
      }
      if (isEvidenceStronger(lbHighestEvidence,comp.getLbIdentCode())){
        lbHighestEvidence = comp.getLbIdentCode();
        lbHighestHit = i;
      }
      if (checkMS2 && comp.getMs2Evidence().size()>0){
        for (LdaLBLASTCompareVO comp2 : comp.getMs2Evidence()){
          int ldaHighestEvidenceMS2 = -4;
          int ldaHighestHitMS2 = 0;
          int lbHighestEvidenceMS2 = -4;
          int lbHighestHitMS2 = 0;
          String key = comp2.getCorrectName();
          if (key.endsWith("_FP")) key = key.substring(0,key.length()-"_FP".length());
          if (key.endsWith("_noMS2")) key = key.substring(0,key.length()-"_noMS2".length());
          if (ldaMS2Ev.containsKey(key)){
            ldaHighestEvidenceMS2 = ldaMS2Ev.get(key);
            ldaHighestHitMS2 = ldaMS2Hit.get(key);
          }
          if (lbMS2Ev.containsKey(key)){
            lbHighestEvidenceMS2 = lbMS2Ev.get(key);
            lbHighestHitMS2 = lbMS2Hit.get(key);
          }
          if (isEvidenceStronger(ldaHighestEvidenceMS2,comp2.getLdaIdentCode())){
            ldaHighestEvidenceMS2 = comp2.getLdaIdentCode();
            ldaHighestHitMS2 = i;
          }
          if (isEvidenceStronger(lbHighestEvidenceMS2,comp2.getLbIdentCode())){
            lbHighestEvidenceMS2 = comp2.getLbIdentCode();
            lbHighestHitMS2 = i;
          }
          ldaMS2Ev.put(key, ldaHighestEvidenceMS2);
          ldaMS2Hit.put(key, ldaHighestHitMS2);
          lbMS2Ev.put(key, lbHighestEvidenceMS2);
          lbMS2Hit.put(key, lbHighestHitMS2);          
        }
      }
    }
    for (int i=0;i!=compares.size();i++){
      LdaLBLASTCompareVO comp = compares.get(i);
      if (i!=ldaHighestHit) comp.setIgnoreLDA(true);
      if (i!=lbHighestHit) comp.setIgnoreLB(true);
      if (checkMS2 && comp.getMs2Evidence().size()>0){
        for (LdaLBLASTCompareVO comp2 : comp.getMs2Evidence()){
          String key = comp2.getCorrectName();
          if (key.endsWith("_FP")) key = key.substring(0,key.length()-"_FP".length());
          if (key.endsWith("_noMS2")) key = key.substring(0,key.length()-"_noMS2".length());
          if (i!=ldaMS2Hit.get(key)) comp2.setIgnoreLDA(true);
          if (i!=lbMS2Hit.get(key)) comp2.setIgnoreLB(true);
        }
      }
    }
  }
  
  private boolean isEvidenceStronger(int current, int ev){
    if (current==-4) return true;
    else if (current==LdaLBLASTCompareVO.NOT_FOUND_OR_MS1_ONLY){ 
      if (ev!=current) return true;
    }else if (current==LdaLBLASTCompareVO.FALSE_NEGATIVE){
      if (ev!=LdaLBLASTCompareVO.NOT_FOUND_OR_MS1_ONLY && ev!=current) return true;
    } else {
      if (ev>current && ev!=LdaLBLASTCompareVO.NOT_FOUND_OR_MS1_ONLY && ev!=LdaLBLASTCompareVO.FALSE_NEGATIVE) return true;
    }
    return false;
  }
  
  private void detectLBNotDetected() throws LipidCombinameEncodingException{
    Hashtable<String,Vector<String>> detectedByLB = new Hashtable<String,Vector<String>>();
    Hashtable<String,Hashtable<String,Vector<SimpleValueObject>>> foundByLDAOnly = new Hashtable<String,Hashtable<String,Vector<SimpleValueObject>>>();
    String baseDir = "L:\\Biological_Experiment_LDA2\\LipidBlast\\";
    Vector<String> excels = new Vector<String>();
    excels.add(baseDir+"positive\\002_liver2-1_Orbitrap_CID_pos_LB10_comp.xlsx");
    excels.add(baseDir+"positive\\003_liver2-1_Orbitrap_CID_pos_LB10_comp.xlsx");
    excels.add(baseDir+"positive\\004_liver2-1_Orbitrap_CID_pos_LB10_comp.xlsx");
    excels.add(baseDir+"positive\\005_liver2-1_Orbitrap_CID_pos_LB10_comp.xlsx");
    excels.add(baseDir+"positive\\006_liver2-1_Orbitrap_CID_pos_LB10_comp.xlsx");

    excels.add(baseDir+"negative\\002_liver2-1_Orbitrap_CID_neg_LB10_comp.xlsx");
    excels.add(baseDir+"negative\\003_liver2-1_Orbitrap_CID_neg_LB10_comp.xlsx");
    excels.add(baseDir+"negative\\004_liver2-1_Orbitrap_CID_neg_LB10_comp.xlsx");
    excels.add(baseDir+"negative\\005_liver2-1_Orbitrap_CID_neg_LB10_comp.xlsx");
    excels.add(baseDir+"negative\\006_liver2-1_Orbitrap_CID_neg_LB10_comp.xlsx");

    for (String excel : excels){
      Vector<Hashtable<String,Vector<String>>> detectedAndNotDetected = detectLBNotDetectedFromExcel(excel);
      Hashtable<String,Vector<String>> ldaOnlyDetected = detectedAndNotDetected.get(0);
      Hashtable<String,Vector<String>> lbDetected = detectedAndNotDetected.get(1);
      Hashtable<String,Vector<String>> ldaOriginal = detectedAndNotDetected.get(2);
      //remove LDA detected and fill up the lbDetected
      for (String className : lbDetected.keySet()){
        Vector<String> lbClassDetected = lbDetected.get(className);
        Hashtable<String,Vector<SimpleValueObject>> ldaClassFound = new Hashtable<String,Vector<SimpleValueObject>>();
        if (foundByLDAOnly.containsKey(className)) ldaClassFound = foundByLDAOnly.get(className);
        Vector<String> alreadyDetected = new Vector<String>();
        if (detectedByLB.containsKey(className)) alreadyDetected = detectedByLB.get(className);
        for (String detected : lbClassDetected){
          Vector<String> keys = new Vector<String>(ldaClassFound.keySet());
          int pos = this.isFACombiInList(detected, keys);
          if (pos>-1) ldaClassFound.remove(keys.get(pos));
          pos = this.isFACombiInList(detected, alreadyDetected);
          if (pos==-1) alreadyDetected.add(detected);
        }
        detectedByLB.put(className, alreadyDetected);
        foundByLDAOnly.put(className, ldaClassFound);
      }
      //add the ones that were found by LDA only
      for (String className : ldaOnlyDetected.keySet()){
        Vector<String> ldaOnlys = ldaOnlyDetected.get(className);
        Vector<String> ldaOriginals = ldaOriginal.get(className);
        
        Vector<String> lbClass = new Vector<String>();
        if (detectedByLB.containsKey(className)) lbClass = detectedByLB.get(className);
        Hashtable<String,Vector<SimpleValueObject>> ldaClassFound = new Hashtable<String,Vector<SimpleValueObject>>();
        if (foundByLDAOnly.containsKey(className)) ldaClassFound = foundByLDAOnly.get(className);
        for (int i=0;i!=ldaOnlys.size();i++){
          String structure = ldaOnlys.get(i);
          int pos = this.isFACombiInList(structure, lbClass);
          if (pos>-1) continue;
          Vector<SimpleValueObject> info = new Vector<SimpleValueObject>();
          Vector<String> keys = new Vector<String>(ldaClassFound.keySet());
          pos = this.isFACombiInList(structure, keys);
          if (pos>-1){
            String key = keys.get(pos);
            info = ldaClassFound.get(key);
            structure = key;
          }
          String fileName = excel.substring(excel.lastIndexOf("\\")+1);
          SimpleValueObject svo = new SimpleValueObject(fileName,ldaOriginals.get(i));
          info.add(svo);
          ldaClassFound.put(structure, info);
        }
        if (ldaClassFound.size()>0) foundByLDAOnly.put(className, ldaClassFound);
      }
    }
    try{
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(baseDir+"LDAFoundOnly.xlsx"));
      Workbook resultWorkbook = new XSSFWorkbook();
      Sheet sheet = resultWorkbook.createSheet("LDA_only");
      int rowCount = 0;
      int sumColumn = 0;
      int structureColumn = 1;
      int detailsColumn = 2;
      for (String className : foundByLDAOnly.keySet()){
        if (rowCount!=0) rowCount++;
        Row row = sheet.createRow(rowCount);
        rowCount++;
        Cell cell = row.createCell(sumColumn);
        cell.setCellType(XSSFCell.CELL_TYPE_STRING);
        cell.setCellValue(className);
        Hashtable<String,Vector<SimpleValueObject>> ldaClassOnly = foundByLDAOnly.get(className);
        Hashtable<Integer,Hashtable<Integer,Integer>> sumSpecies = new Hashtable<Integer,Hashtable<Integer,Integer>>();
        for (String structure : ldaClassOnly.keySet()){
          int[] cAndDbs = getCAtomsAndDoubleBonds(structure);
          Hashtable<Integer,Integer> dbs = new Hashtable<Integer,Integer>();
          if (sumSpecies.containsKey(cAndDbs[0])) dbs = sumSpecies.get(cAndDbs[0]);
          dbs.put(cAndDbs[1], cAndDbs[1]);
          sumSpecies.put(cAndDbs[0], dbs);
        }
        List<Integer> cAtoms = new  ArrayList<Integer>(sumSpecies.keySet());
        Collections.sort(cAtoms);
        for (Integer cAtom : cAtoms){
          List<Integer> dbs = new ArrayList<Integer>(sumSpecies.get(cAtom).keySet());
          Collections.sort(dbs);
          for (Integer db : dbs){
            boolean firstSumStructure = true;
            for (String structure : ldaClassOnly.keySet()){
              int[] cAndDbs = getCAtomsAndDoubleBonds(structure);
              if (cAtom!=cAndDbs[0] || db!=cAndDbs[1]) continue;
              row = sheet.createRow(rowCount);
              rowCount++;
              if (firstSumStructure){
                cell = row.createCell(sumColumn);
                cell.setCellType(XSSFCell.CELL_TYPE_STRING);
                cell.setCellValue(className+" "+cAtom+":"+db);                
              }
              Vector<SimpleValueObject> vos = ldaClassOnly.get(structure);
              Hashtable<String,Integer> mostOftenOccuringStructure = new Hashtable<String,Integer>();
              String detailsString = "";
              for (SimpleValueObject vo : vos){
                if (detailsString.length()>0) detailsString+="; ";
                detailsString += vo.getLabel().substring(0,4);
                if (vo.getLabel().indexOf("pos")>-1) detailsString += "pos";
                else if (vo.getLabel().indexOf("neg")>-1) detailsString += "neg";
                detailsString+="="+vo.getValue();
                int count = 0;
                if (mostOftenOccuringStructure.containsKey(vo.getValue())) count = mostOftenOccuringStructure.get(vo.getValue());
                count++;
                mostOftenOccuringStructure.put(vo.getValue(),count);
              }
              String mostPresentStructure = structure;
              int mostOften = 0;
              for (String struct : mostOftenOccuringStructure.keySet()){
                if (mostOftenOccuringStructure.get(struct)>mostOften){
                  mostOften = mostOftenOccuringStructure.get(struct);
                  mostPresentStructure = struct;
                }
              }              
              firstSumStructure = false;
              cell = row.createCell(structureColumn);
              cell.setCellType(XSSFCell.CELL_TYPE_STRING);
              cell.setCellValue(mostPresentStructure);
              cell = row.createCell(detailsColumn);
              cell.setCellType(XSSFCell.CELL_TYPE_STRING);
              cell.setCellValue(detailsString);
            }    
          }
        }
      }
      resultWorkbook.write(out);
      out.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }
  }
  
  private Vector<Hashtable<String,Vector<String>>> detectLBNotDetectedFromExcel(String excel){
    Vector<Hashtable<String,Vector<String>>> detectedAndNotDetected = new Vector<Hashtable<String,Vector<String>>>();
    Hashtable<String,Vector<String>> lbDetected = new Hashtable<String,Vector<String>>();
    Hashtable<String,Vector<String>> ldaOnlyDetected = new Hashtable<String,Vector<String>>();
    Hashtable<String,Vector<String>> ldaOriginal = new Hashtable<String,Vector<String>>();
    try{
      InputStream myxls = new FileInputStream(excel);
      Workbook workbook = new XSSFWorkbook(excel);
      for (int i=0; i!=workbook.getNumberOfSheets();i++){
        Sheet sheet = workbook.getSheetAt(i);
        String className = sheet.getSheetName();
        //if (!className.equalsIgnoreCase("TG")) continue;
        Row row = sheet.getRow(0);
        int structureColumn = -1;
        int ldaStructureColumn = -1;
        int ldaCodeColumn = -1;
        int lbCodeColumn = -1;
        int nameColumn = -1;
        for (int j=0; j!=row.getLastCellNum(); j++){
          Cell cell = row.getCell(j);
          String contents = "";
          int cellType = -1;
          if (cell==null) continue;
          cellType = cell.getCellType();
          if (cellType!=XSSFCell.CELL_TYPE_STRING) continue;
          contents = cell.getStringCellValue();
          if (contents.equalsIgnoreCase("Structure")) structureColumn = j;
          else if (contents.equalsIgnoreCase("LDA")) ldaStructureColumn = j;
          else if (contents.equalsIgnoreCase("LDA-Code")) ldaCodeColumn = j;
          else if (contents.equalsIgnoreCase("LB-Code")) lbCodeColumn = j;
          else if (contents.equalsIgnoreCase("Name")) nameColumn = j;
        }
        if (ldaCodeColumn>-1 && lbCodeColumn>-1 && structureColumn<0) structureColumn = nameColumn;
        if (ldaCodeColumn<0 && lbCodeColumn<0 && structureColumn<0) continue;
        
        Vector<String> lbClassDetected = new Vector<String>();
        Vector<String> ldaClassDetected = new Vector<String>();
        Vector<String> ldaOriginalClassDetected = new Vector<String>();
        for (int rowNum=1; rowNum<sheet.getLastRowNum()+1;rowNum++){
          row = sheet.getRow(rowNum);
          String structureOriginal = "";
          String structure = "";
          String ldaStructure = "";
          int ldaCode = -1000;
          int lbCode = -1000;
          
          if (row==null) continue;
          //reading structureInformation
          Cell cell = row.getCell(structureColumn);
          if (cell==null || cell.getCellType()!=XSSFCell.CELL_TYPE_STRING) continue;
          structureOriginal = cell.getStringCellValue().trim();
          if (structureOriginal.endsWith("_FP") || structureOriginal.endsWith("_noMS2") || structureOriginal.length()==0) continue;
          structure = structureOriginal.replaceAll("/", "_");
          
          cell = row.getCell(ldaStructureColumn);
          if (cell==null || cell.getCellType()!=XSSFCell.CELL_TYPE_STRING) continue;
          ldaStructure = cell.getStringCellValue().trim();
          
         
          //reading ldaCode
          cell = row.getCell(ldaCodeColumn);
          if (cell==null) continue;
          int cellType = cell.getCellType();
          if (cellType==XSSFCell.CELL_TYPE_STRING)
            try{ldaCode = Integer.parseInt(cell.getStringCellValue());}catch(NumberFormatException nfx){continue;}
          else if (cellType==XSSFCell.CELL_TYPE_NUMERIC || cellType==XSSFCell.CELL_TYPE_FORMULA){
            try{ldaCode = (new Double(cell.getNumericCellValue()).intValue());}catch(NumberFormatException nfx){continue;}
          } else continue;
          
          //reading ldaCode
          cell = row.getCell(lbCodeColumn);
          if (cell==null) continue;
          cellType = cell.getCellType();
          if (cellType==XSSFCell.CELL_TYPE_STRING)
            try{lbCode = Integer.parseInt(cell.getStringCellValue());}catch(NumberFormatException nfx){continue;}
          else if (cellType==XSSFCell.CELL_TYPE_NUMERIC || cellType==XSSFCell.CELL_TYPE_FORMULA){
            try{lbCode = (new Double(cell.getNumericCellValue()).intValue());}catch(NumberFormatException nfx){continue;}
          } else continue;
          
          if (ldaCode==0 && lbCode==0) continue;
          
          if (lbCode>0){
            int pos = isFACombiInList(structure, ldaClassDetected);
            if (pos>-1){
              ldaClassDetected.remove(pos);
              ldaOriginalClassDetected.remove(pos);
            }
            pos = isFACombiInList(structure, lbClassDetected);
            if (pos==-1) lbClassDetected.add(structure);
          } else if (ldaCode>0 && (lbCode==-1 || lbCode==-3)){
            int pos = isFACombiInList(structure, lbClassDetected);
            if (pos>-1) continue;
            pos = isFACombiInList(structure, ldaClassDetected);
            if (pos>-1) continue;
            ldaClassDetected.add(structure);
            ldaOriginalClassDetected.add(ldaStructure);
          }
          //System.out.println(className+" "+structureOriginal+";"+ldaCode+";"+lbCode);
        }
        if (lbClassDetected.size()>0 || ldaClassDetected.size()>0){
          lbDetected.put(className, lbClassDetected);
          ldaOnlyDetected.put(className, ldaClassDetected);
          ldaOriginal.put(className, ldaOriginalClassDetected);
        }
      }

    } catch (Exception ex){
      ex.printStackTrace();
    }
    detectedAndNotDetected.add(ldaOnlyDetected);
    detectedAndNotDetected.add(lbDetected);
    detectedAndNotDetected.add(ldaOriginal);
    return detectedAndNotDetected;
  }
  
  private int isFACombiInList(String structure, Vector<String> list) throws LipidCombinameEncodingException{
    for (int i=0; i!=list.size(); i++){
      if (StaticUtils.isAPermutedVersion(structure, list.get(i),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) return i;
    }
    return -1;
  }
  
  private int[] getCAtomsAndDoubleBonds(String structure){
    int[] cAndDbs = new int[2];
    String[] fas = structure.split("_");
    cAndDbs[0] = 0;
    cAndDbs[1] = 0;
    for (String fa:fas){
      if (fa.equalsIgnoreCase("-")) continue;
      String toSplit = new String(fa);
      if (fa.startsWith("P-")) toSplit = toSplit.substring(2);
      cAndDbs[0] += Integer.parseInt(toSplit.substring(0,toSplit.indexOf(":")));
      cAndDbs[1] += Integer.parseInt(toSplit.substring(toSplit.indexOf(":")+1));
    }
    return cAndDbs;
  }
  
  private void tryXmlStax(){
    long time = System.currentTimeMillis();
    /*    AddScan[] adders = new AddScan[1];
    adders[0] = this;
    CgReader reader = new MzXmlReader(adders,true,1000);
    try {
      reader.ReadFile("D:\\Evelyn\\20171006\\SRM1950_1.mzXML");
      //reader.ReadFile("D:\\Kristaps\\20171129\\TG quant NIST\\MCC007_Lipid01_NIST1_20171124.mzXML",true);
      
      System.out.println(reader.getLowestMz()+";"+reader.getHighestMz());
      System.out.println(reader.usesPolaritySwitching());
    }
    catch (Exception e) {
      e.printStackTrace();
    }*/
    
    //String filePath = "D:\\testMzXML\\bruker\\Mix25_RD1_01_1360.mzXML";
    String filePath = "D:\\testMzXML\\thermo\\Schev_Meth_test_mono-iso-CID_120731120418.mzXML";
    //String filePath = "D:\\testMzXML\\waters\\20151103_GNR_LDA_Exp3_1_01.mzXML";
    //String filePath = "D:\\testMzXML\\absciex\\Data20141110_versch_CE_pos+neg1-022b.mzXML";
    //String filePath = "D:\\testMzXML\\agilent\\Sample1_pos_01.mzXML";
//    String filePath = "D:\\positionIsomers\\test\\Mix_5uM-95min-1.mzXML";
    //String filePath = "D:\\Kristaps\\20171129\\TG quant NIST\\MCC007_Lipid01_NIST1_20171124.mzXML";
    //System.out.println(this.getPiecesForChromTranslation(filePath));
//    RawToChromTranslator translator = new RawToChromTranslator(filePath,"mzXML",600,7,
//        100,1,false);
    RawToChromTranslator translator = new RawToChromTranslator(filePath,"mzXML",1,7,
        1000,1,true);
//    RawToChromTranslator translator = new RawToChromTranslator(filePath,"mzXML",10,7,
//        1000,1,true);
//    RawToChromTranslator translator = new RawToChromTranslator(filePath,"mzXML",50,7,
//        1000,5,true);
//    RawToChromTranslator translator = new RawToChromTranslator(filePath,"mzXML",50,7,
//        1000,1,true);

//    String filePath = "D:\\Evelyn\\20171006\\SRM1950_1.mzXML";
//    RawToChromTranslator translator = new RawToChromTranslator(filePath,"mzXML",100,7,
//        1000,1,true);

//    RawToChromTranslator translator = new RawToChromTranslator(filePath,"mzXML",800,7,
//        1000,1,true);
    
    try {
      translator.translateToChromatograms();
      Set<String> fileNames = translator.getOriginalFileNames();
//      for (String name :fileNames){
//        System.out.println("fileName: "+name);
//      }
    }
    catch (CgException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
    /****XMLInputFactory factory = XMLInputFactory.newInstance();
    factory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, true);
    XMLStreamReader streamReader;
    try {
      streamReader = factory.createXMLStreamReader(new FileReader("D:\\testMzXML\\bruker\\Mix25_RD1_01_1360.mzXML"));
      while(streamReader.hasNext()){
        int eventType = streamReader.next();

        if(eventType == XMLStreamReader.START_ELEMENT){
           // System.out.println(streamReader.getLocalName());
        }
    }

      
      streamReader.close();
    }
    catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch (XMLStreamException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }*/
    
    System.out.println((System.currentTimeMillis()-time)/1000l+" sec");
  }

  @Override
  public void AddHeader(CgScanHeader arg0) throws CgException
  {
    // TODO Auto-generated method stub
    
  }

  @Override
  public void AddScan(CgScan arg0) throws CgException
  {
    // TODO Auto-generated method stub
    
  }

  @Override
  public void addParentFileName(String arg0) throws CgException
  {
    // TODO Auto-generated method stub
    
  }

  @Override
  public CgScan getLastBaseScan()
  {
    // TODO Auto-generated method stub
    return null;
  }

  @Override
  public void setStartStopHeader(CgScanHeader arg0) throws CgException
  {
    // TODO Auto-generated method stub
    
  }
  
  
  private void getValidOrbitrapCIDSpeciesPositive(LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses,
      Hashtable<String,LipidClassInfoVO> lipidClassInfo, LinkedHashMap<String,Boolean> adducts){
    lipidClasses.put("P-PC", FoundBiologicalSpecies.getPPCSpeciesOrbitrap());
    adducts.put("H", true);
    adducts.put("Na", true);
    lipidClassInfo.put("P-PC", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("P-PE", FoundBiologicalSpecies.getPPESpeciesOrbitrap());
    adducts.put("H", true);
    adducts.put("Na", true);
    lipidClassInfo.put("P-PE", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("LPE", FoundBiologicalSpecies.getLPESpeciesOrbitrap());
    adducts.put("H", false);
    adducts.put("Na", false);
    lipidClassInfo.put("LPE", new LipidClassInfoVO(1,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PS", FoundBiologicalSpecies.getPSSpeciesOrbitrap());
    adducts.put("H", true);
    lipidClassInfo.put("PS", new LipidClassInfoVO(2,false,-1,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PC", FoundBiologicalSpecies.getPCSpeciesOrbitrap());
    adducts.put("H", true);
    adducts.put("Na", true);
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClassInfo.put("PC", new LipidClassInfoVO(2,true,0.7d,adducts));
    lipidClasses.put("PE", FoundBiologicalSpecies.getPESpeciesOrbitrap());
    adducts.put("H", true);
    adducts.put("Na", true);
    lipidClassInfo.put("PE", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("Cer", FoundBiologicalSpecies.getCerSpeciesOrbitrap());
    adducts.put("H", false);
    lipidClassInfo.put("Cer", new LipidClassInfoVO(1,true,0.7d,adducts));
    
    //this has to be before LPC, DG, TG and SM
    correctRetentionTimes(-0.2d,lipidClasses);
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("LPC", FoundBiologicalSpecies.getLPCSpeciesOrbitrap());
    adducts.put("H", false);
    adducts.put("Na", false);
    lipidClassInfo.put("LPC", new LipidClassInfoVO(1,true,1.2d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("DG", FoundBiologicalSpecies.getDGSpeciesOrbitrap());
    adducts.put("Na", true);
    adducts.put("NH4", true);
    lipidClassInfo.put("DG", new LipidClassInfoVO(3,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("TG", FoundBiologicalSpecies.getTGSpeciesOrbitrap());
    adducts.put("NH4", true);
    adducts.put("Na", true);
    lipidClassInfo.put("TG", new LipidClassInfoVO(3,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("SM", FoundBiologicalSpecies.getSMSpeciesOrbitrap());
    adducts.put("H", false);
    adducts.put("Na", false);
    lipidClassInfo.put("SM", new LipidClassInfoVO(1,true,0.7d,adducts));
  }
  
  private void getValid4000QTRAPSpeciesPositive(LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses,
      Hashtable<String,LipidClassInfoVO> lipidClassInfo, LinkedHashMap<String,Boolean> adducts){
    //it is not recommended to analyzer P-PC with LDA on a 4000 QTRAP
//    lipidClasses.put("P-PC", FoundBiologicalSpecies.getPPCSpecies4000QTRAP());
//    adducts.put("H", "H");
//    adducts.put("Na", "Na");
//    lipidClassInfo.put("P-PC", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("P-PE", FoundBiologicalSpecies.getPPESpecies4000QTRAP());
    adducts.put("H", true);
    adducts.put("Na", true);
    lipidClassInfo.put("P-PE", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("LPE", FoundBiologicalSpecies.getLPESpecies4000QTRAP());
    adducts.put("H", false);
    adducts.put("Na", false);
    lipidClassInfo.put("LPE", new LipidClassInfoVO(1,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PS", FoundBiologicalSpecies.getPSSpecies4000QTRAP());
    adducts.put("H", true);
    lipidClassInfo.put("PS", new LipidClassInfoVO(2,false,-1,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PC", FoundBiologicalSpecies.getPCSpecies4000QTRAP());
    adducts.put("H", true);
    adducts.put("Na", true);
    lipidClassInfo.put("PC", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PE", FoundBiologicalSpecies.getPESpecies4000QTRAP());
    adducts.put("H", true);
    adducts.put("Na", true);
    lipidClassInfo.put("PE", new LipidClassInfoVO(2,true,0.7d,adducts));
    
    //this has to be before LPC, DG, TG and SM
    ////correctRetentionTimes(-0.2d,lipidClasses);
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("LPC", FoundBiologicalSpecies.getLPCSpecies4000QTRAP());
    adducts.put("H", false);
    adducts.put("Na", false);
    lipidClassInfo.put("LPC", new LipidClassInfoVO(1,true,1.2d,adducts));
    //it is not recommended to analyze DG with a 4000 QTRAP
//    adducts = new LinkedHashMap<String,String>();
//    lipidClasses.put("DG", FoundBiologicalSpecies.getDGSpecies4000QTRAP());
//    adducts.put("Na", "Na");
//    lipidClassInfo.put("DG", new LipidClassInfoVO(3,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("TG", FoundBiologicalSpecies.getTGSpecies4000QTRAP());
    adducts.put("NH4", true);
    adducts.put("Na", true);
    lipidClassInfo.put("TG", new LipidClassInfoVO(3,true,0.6d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("SM", FoundBiologicalSpecies.getSMSpecies4000QTRAP());
    adducts.put("H", true);
    adducts.put("Na", true);
    lipidClassInfo.put("SM", new LipidClassInfoVO(1,true,0.7d,adducts));

  }

  private void getValidOrbitrapCIDSpeciesNegative(LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses,
      Hashtable<String,LipidClassInfoVO> lipidClassInfo, LinkedHashMap<String,Boolean> adducts){
    lipidClasses.put("PI", FoundBiologicalSpecies.getPISpeciesOrbitrap());
    adducts.put("-H", true);
    lipidClassInfo.put("PI", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("P-PE", FoundBiologicalSpecies.getPPESpeciesOrbitrap());
    adducts.put("-H", true);
    lipidClassInfo.put("P-PE", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("LPE", FoundBiologicalSpecies.getLPESpeciesOrbitrap());
    adducts.put("-H", false);
    lipidClassInfo.put("LPE", new LipidClassInfoVO(1,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PS", FoundBiologicalSpecies.getPSSpeciesOrbitrap());
    adducts.put("-H", true);
    lipidClassInfo.put("PS", new LipidClassInfoVO(2,false,-1,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PC", FoundBiologicalSpecies.getPCSpeciesOrbitrap());
    adducts.put("HCOO", true);
    adducts.put("-CH3", true);
    lipidClassInfo.put("PC", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PE", FoundBiologicalSpecies.getPESpeciesOrbitrap());
    adducts.put("-H", true);
    lipidClassInfo.put("PE", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PG", FoundBiologicalSpecies.getPGSpeciesOrbitrap());
    adducts.put("-H", true);
    lipidClassInfo.put("PG", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("Cer", FoundBiologicalSpecies.getCerSpeciesOrbitrap());
    adducts.put("-H", false);
    lipidClassInfo.put("Cer", new LipidClassInfoVO(1,true,0.7d,adducts));
  }
  
  private void getValid4000QTRAPSpeciesNegative(LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses,
      Hashtable<String,LipidClassInfoVO> lipidClassInfo, LinkedHashMap<String,Boolean> adducts){
    lipidClasses.put("PI", FoundBiologicalSpecies.getPISpecies4000QTRAP());
    adducts.put("-H", true);
    lipidClassInfo.put("PI", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("P-PE", FoundBiologicalSpecies.getPPESpecies4000QTRAP());
    adducts.put("-H", true);
    lipidClassInfo.put("P-PE", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("LPE", FoundBiologicalSpecies.getLPESpecies4000QTRAP());
    adducts.put("-H", false);
    lipidClassInfo.put("LPE", new LipidClassInfoVO(1,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PS", FoundBiologicalSpecies.getPSSpecies4000QTRAP());
    adducts.put("-H", true);
    lipidClassInfo.put("PS", new LipidClassInfoVO(2,true,0.7,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PC", FoundBiologicalSpecies.getPCSpecies4000QTRAP());
    adducts.put("HCOO", true);
    adducts.put("-CH3", true);
    lipidClassInfo.put("PC", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PE", FoundBiologicalSpecies.getPESpecies4000QTRAP());
    adducts.put("-H", true);
    lipidClassInfo.put("PE", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("PG", FoundBiologicalSpecies.getPGSpecies4000QTRAP());
    adducts.put("-H", true);
    lipidClassInfo.put("PG", new LipidClassInfoVO(2,true,0.7d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("Cer", FoundBiologicalSpecies.getCerSpecies4000QTRAP());
    adducts.put("-H", false);
    lipidClassInfo.put("Cer", new LipidClassInfoVO(1,true,0.7d,adducts));
  }
  
  private void generateDetailsBiologicalExperiment(){
    Vector<String> quantitationFiles = new Vector<String>();
    quantitationFiles.add("D:\\BiologicalExperiment\\massLists\\positive\\positive.xlsx");
    quantitationFiles.add("D:\\BiologicalExperiment\\massLists\\negative\\negative.xlsx");

    // these are the values for the G6550A QTOF
//    String positiveListFile = "D:\\BiologicalExperiment\\Singapore\\SpeciesDetectable_MSMS_Singapore.xlsx";
//    String resultFile = "D:\\BiologicalExperiment\\Singapore\\SpeciesDetectable_MSMS_G6550A_QTOF_Details_generated.xlsx";
//    String sheetName = "G6550A QTOF";
//        
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\positive\\Sample2.1_pos_01_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\positive\\Sample2.1_pos_02_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\positive\\Sample2.1_pos_03_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\positive\\Sample2.1_pos_04_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\positive\\Sample2.1_pos_05_positive.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//    
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\negative\\Sample2.1_neg_01_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\negative\\Sample2.1_neg_02_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\negative\\Sample2.1_neg_03_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\negative\\Sample2.1_neg_04_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Singapore\\negative\\Sample2.1_neg_05_negative.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//    
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.004d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getQTOFG6550Modifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getQTOFG6550PreElutionRanges();

    // these are the values for the Orbitrap Elite
//    String positiveListFile = "D:\\BiologicalExperiment\\Cambridge\\SpeciesDetectable_MSMS_Cambridge.xlsx";
//    String resultFile = "D:\\BiologicalExperiment\\Cambridge\\SpeciesDetectable_MSMS_Orbitrap_Elite_Details_generated.xlsx";
//    String sheetName = "Orbitrap Elite";
//        
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\positive\\3a_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\positive\\3b_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\positive\\3c_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\positive\\3d_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\positive\\3e_positive.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//    
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\negative\\3a_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\negative\\3b_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\negative\\3c_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\negative\\3d_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Cambridge\\negative\\3e_negative.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//    
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.003d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapEliteModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getOrbitrapElitePreElutionRanges();

    
    //these are the values for the Orbitrap_CID
//    String positiveListFile = "D:\\BiologicalExperiment\\Orbitrap_CID\\SpeciesDetectable_MSMS_Orbitrap_CID.xlsx";
//    String resultFile = "D:\\BiologicalExperiment\\Orbitrap_CID\\SpeciesDetectable_MSMS_Orbitrap_CID_Details_generated.xlsx";
//    String sheetName = "Orbitrap Velos Pro CID";
//    
//    
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\positive\\002_liver2-1_Orbitrap_CID_pos_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\positive\\003_liver2-1_Orbitrap_CID_pos_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\positive\\004_liver2-1_Orbitrap_CID_pos_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\positive\\005_liver2-1_Orbitrap_CID_pos_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\positive\\006_liver2-1_Orbitrap_CID_pos_positive.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//    
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\negative\\002_liver2-1_Orbitrap_CID_neg_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\negative\\003_liver2-1_Orbitrap_CID_neg_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\negative\\004_liver2-1_Orbitrap_CID_neg_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\negative\\005_liver2-1_Orbitrap_CID_neg_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_CID\\negative\\006_liver2-1_Orbitrap_CID_neg_negative.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//      
//    
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.002d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapCIDModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getOrbitrapCIDPreElutionRanges();
    
    // these are the values for the Orbitrap_HCD
    String positiveListFile = "D:\\BiologicalExperiment\\Orbitrap_HCD\\SpeciesDetectable_MSMS_Orbitrap_HCD.xlsx";
    String resultFile = "D:\\BiologicalExperiment\\Orbitrap_HCD\\SpeciesDetectable_MSMS_Orbitrap_HCD_Details_generated.xlsx";
    String sheetName = "Orbitrap Velos Pro HCD";
      
    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
    Vector<String> filesOfOneIonMode = new Vector<String>();
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\positive\\020_liver2-1_Orbitrap_HCD_pos_positive.xlsx");
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\positive\\021_liver2-1_Orbitrap_HCD_pos_positive.xlsx");
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\positive\\022_liver2-1_Orbitrap_HCD_pos_positive.xlsx");
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\positive\\023_liver2-1_Orbitrap_HCD_pos_positive.xlsx");
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\positive\\024_liver2-1_Orbitrap_HCD_pos_positive.xlsx");
    ldaFilesIonMode.put("+", filesOfOneIonMode);
  
    filesOfOneIonMode = new Vector<String>();
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\negative\\020_liver2-1_Orbitrap_HCD_neg_negative.xlsx");
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\negative\\021_liver2-1_Orbitrap_HCD_neg_negative.xlsx");
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\negative\\022_liver2-1_Orbitrap_HCD_neg_negative.xlsx");
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\negative\\023_liver2-1_Orbitrap_HCD_neg_negative.xlsx");
    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Orbitrap_HCD\\negative\\024_liver2-1_Orbitrap_HCD_neg_negative.xlsx");
    ldaFilesIonMode.put("-", filesOfOneIonMode);
  
    int decimalPlacesForMz = 4;
    double highestMzStdevAllowed = 0.002d;
    double highestRtStdevAllowed = 0.2;
    double highestAbundanceCVAllowed = 50;
    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapHCDModifications();
    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getOrbitrapHCDPreElutionRanges();


    // these are the values for the Q Exactive
//    String positiveListFile = "D:\\BiologicalExperiment\\Stockholm\\SpeciesDetectable_MSMS_Stockholm.xlsx";
//    String resultFile = "D:\\BiologicalExperiment\\Stockholm\\SpeciesDetectable_MSMS_QExactive_Details_generated.xlsx";
//    String sheetName = "Q Exactive";
//      
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\positive\\20150709_Liver1_pos-01_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\positive\\20150709_Liver1_pos-02_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\positive\\20150709_Liver1_pos-03_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\positive\\20150709_Liver1_pos-04_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\positive\\20150709_Liver1_pos-05_positive.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//  
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\negative\\20150709_Liver1_neg-01_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\negative\\20150709_Liver1_neg-02_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\negative\\20150709_Liver1_neg-03_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\negative\\20150709_Liver1_neg-04_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\Stockholm\\negative\\20150709_Liver1_neg-05_negative.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//  
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.006d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getQExactiveModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getQExactivePreElutionRanges();

    // these are the values for the 4000 QTRAP
//    String positiveListFile = "D:\\BiologicalExperiment\\QTRAP\\SpeciesDetectable_MSMS_QTRAP.xlsx";
//    String resultFile = "D:\\BiologicalExperiment\\QTRAP\\SpeciesDetectable_MSMS_4000QTRAP_Details_generated.xlsx";
//    String sheetName = "4000 QTRAP";
//      
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\positive\\Data20151002_QTrap_Liver-002_QTrap_Liver1-1_pos_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\positive\\Data20151002_QTrap_Liver-003_QTrap_Liver1-1_pos_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\positive\\Data20151002_QTrap_Liver-004_QTrap_Liver1-1_pos_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\positive\\Data20151002_QTrap_Liver-005_QTrap_Liver1-1_pos_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\positive\\Data20151002_QTrap_Liver-006_QTrap_Liver1-1_pos_positive.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//  
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\negative\\Data20151002_QTrap_Liver-021_QTrap_Liver1-1_neg_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\negative\\Data20151002_QTrap_Liver-022_QTrap_Liver1-1_neg_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\negative\\Data20151002_QTrap_Liver-023_QTrap_Liver1-1_neg_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\negative\\Data20151002_QTrap_Liver-024_QTrap_Liver1-1_neg_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\QTRAP\\negative\\Data20151002_QTrap_Liver-025_QTrap_Liver1-1_neg_negative.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//  
//    int decimalPlacesForMz = 2;
//    double highestMzStdevAllowed = 0.15d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.get4000QTRAPModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.get4000QTRAPPreElutionRanges();

    // these are the values for the QTRAP 6500
//    String positiveListFile = "D:\\BiologicalExperiment\\SanDiego\\SpeciesDetectable_MSMS_SanDiego.xlsx";
//    String resultFile = "D:\\BiologicalExperiment\\SanDiego\\SpeciesDetectable_MSMS_QTRAP6500_Details_generated.xlsx";
//    String sheetName = "QTRAP 6500";
//      
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\positive\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 pos_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\positive\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 pos (2)_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\positive\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 pos (3)_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\positive\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 pos (4)_positive.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\positive\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 pos (5)_positive.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//  
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\negative\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 neg_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\negative\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 neg (2)_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\negative\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 neg (3)_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\negative\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 neg (4)_negative.xlsx");
//    filesOfOneIonMode.add("D:\\BiologicalExperiment\\SanDiego\\negative\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Liver Extract 1 neg (5)_negative.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//  
//    int decimalPlacesForMz = 2;
//    double highestMzStdevAllowed = 0.15d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getQTRAP6500Modifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getQTRAP6500PreElutionRanges();

    
    
    BufferedOutputStream out = null;
    Hashtable<String,String> appliedMzCorrection = new Hashtable<String,String>();
    int highestNumberOfFiles = 0;
    if (ldaFilesIonMode.get("+").size()>highestNumberOfFiles) highestNumberOfFiles = ldaFilesIonMode.get("+").size();
    if (ldaFilesIonMode.get("-").size()>highestNumberOfFiles) highestNumberOfFiles = ldaFilesIonMode.get("-").size();
    
    try{
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> theoreticalMasses = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();          
      for (String quantFile : quantitationFiles){
        Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> oneFile = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)QuantificationThread.parseQuantExcelFile(quantFile,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f,false).get(3);
        ////Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> oneFile = null;
        for (String className : oneFile.keySet()){
          Hashtable<String,Hashtable<String,QuantVO>> oneClass= oneFile.get(className);
          Hashtable<String,Hashtable<String,QuantVO>> ofClass = new Hashtable<String,Hashtable<String,QuantVO>>();
          if (theoreticalMasses.containsKey(className)) ofClass = theoreticalMasses.get(className);
          for (String analyte:oneClass.keySet()){
            Hashtable<String,QuantVO> oneAnalyte = oneClass.get(analyte);
            Hashtable<String,QuantVO> ofAnalyte = new Hashtable<String,QuantVO>();
            if (ofClass.containsKey(analyte)) ofAnalyte = ofClass.get(analyte);
            for (String mod:oneAnalyte.keySet()) ofAnalyte.put(mod, oneAnalyte.get(mod));
            ofClass.put(analyte, ofAnalyte);
          }
          theoreticalMasses.put(className,ofClass);
        }
      }

      LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> positiveList = this.readPositiveList(positiveListFile);
//      for (String className : positiveList.keySet()){
//        System.out.println(className);
//        for (String species : positiveList.get(className).keySet()){
//          System.out.println("     "+species);
//          for (String structure : positiveList.get(className).get(species).values()){
//            System.out.println("            "+structure);
//          }
//        }
//      }
      //this hashtable is for checking if all of the species are in the positive list
      LinkedHashMap<String,LinkedHashMap<String,Vector<LipidParameterSet>>> checkPositiveList = new LinkedHashMap<String,LinkedHashMap<String,Vector<LipidParameterSet>>>();
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>> results = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>>();
      for (String polarity : ldaFilesIonMode.keySet()){
        for (String ldaFile : ldaFilesIonMode.get(polarity)){
          String fileName = ldaFile.substring(ldaFile.lastIndexOf("\\")+1);
          QuantificationResult quantResult = LDAResultReader.readResultFile(ldaFile,  new Hashtable<String,Boolean>());
          if (quantResult.getConstants().getShift()==0d) appliedMzCorrection.put(fileName, "0");
          else{
            System.out.println(quantResult.getConstants().getShift());
            appliedMzCorrection.put(fileName, Calculator.FormatNumberToString(quantResult.getConstants().getShift(),decimalPlacesForMz));
          }
          Hashtable<String,Vector<LipidParameterSet>> idents = quantResult.getIdentifications();
          //this is for checking if all classes are present
          for (String className : idents.keySet()){
            Vector<LipidParameterSet> classResult = idents.get(className);
            //this is for checking if all classes are present
            boolean found = false;
            for (LipidParameterSet set : classResult){
              if (set instanceof LipidomicsMSnSet){
                found = true;
                break;
              }
            }
            if (found && !positiveList.containsKey(className)){
              System.out.println("The positive list does not contain the lipid class "+className+"! But we identified it in "+fileName);
            }
          // end of checking if all classes are present
          
          //this is for checking if all of the species are in the positive list
//          LinkedHashMap<String,Vector<LipidParameterSet>> classCheck = new LinkedHashMap<String,Vector<LipidParameterSet>>();
//          if (checkPositiveList.containsKey(className)) classCheck = checkPositiveList.get(className);
//          classCheck.put(fileName, classResult);
//          checkPositiveList.put(className, classCheck);
            // end of checking if all of the species are in the positive list

          }
        
          Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> resultsOfExp = new Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>();
          for (String className : positiveList.keySet()){
            if (!idents.containsKey(className)) continue;
            Vector<LipidParameterSet> classResult = idents.get(className);
            LinkedHashMap<String,LinkedHashMap<String,String>> posAnalytes = positiveList.get(className);
            Hashtable<String,Vector<LipidomicsMSnSet>> msnFound = new Hashtable<String,Vector<LipidomicsMSnSet>>();
            for (LipidParameterSet set : classResult){
              if (!(set instanceof LipidomicsMSnSet) || !posAnalytes.containsKey(set.getNameStringWithoutRt())) continue;
              float rt = Float.parseFloat(set.getRt());
              if (sheetName.equalsIgnoreCase("Q Exactive") && (className.equalsIgnoreCase("TG"))){
                if (set.getNameStringWithoutRt().equalsIgnoreCase("58:5") && rt<41.7f) continue;
                else if (set.getNameStringWithoutRt().equalsIgnoreCase("58:7") && rt<40.45f) continue;
                
              }
                  
              Vector<LipidomicsMSnSet> ofOneAnalyte = new Vector<LipidomicsMSnSet>();
              if (msnFound.containsKey(set.getNameStringWithoutRt())) ofOneAnalyte = msnFound.get(set.getNameStringWithoutRt());
                ofOneAnalyte.add((LipidomicsMSnSet)set);
                msnFound.put(set.getNameStringWithoutRt(), ofOneAnalyte);
            }
            resultsOfExp.put(className, msnFound);
          }
          results.put(fileName, resultsOfExp);
        }

      }
      
      
      Workbook resultWorkbook = createDetailsExcelFile(sheetName,highestNumberOfFiles,decimalPlacesForMz,highestMzStdevAllowed,highestRtStdevAllowed,
          highestAbundanceCVAllowed,ldaFilesIonMode,appliedMzCorrection,mods,preElutionRanges,positiveList,theoreticalMasses,results,false,false);
      out = new BufferedOutputStream(new FileOutputStream(resultFile));
      resultWorkbook.write(out);
      resultWorkbook.close();

      // this is for checking the positive list
//      for (String className : checkPositiveList.keySet()){
//        LinkedHashMap<String,Vector<LipidParameterSet>> result = checkPositiveList.get(className);
//        for (String fileName : result.keySet()){
//          Vector<LipidParameterSet> classResults = result.get(fileName);
//          for (LipidParameterSet set : classResults){
//            if (!(set instanceof LipidomicsMSnSet)) continue;
//            LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
//            for (Object msnNames : msn.getMSnIdentificationNames()){    
//              String nameString = "";
//              String faId = "";
//              if (msnNames instanceof Vector){
//                for (String name : (Vector<String>)msnNames){
//                  nameString += name+"|";
//                }
//                faId = nameString.substring(0,nameString.indexOf("|"));
//                nameString = nameString.substring(0,nameString.length()-1);
//              }else{
//                nameString = (String)msnNames;
//                faId = nameString;
//              }
//              LinkedHashMap<String,LinkedHashMap<String,String>> posOfClass = positiveList.get(className);
//              if (msn.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED && !faId.equalsIgnoreCase(msn.getNameStringWithoutRt())){
//                if (!posOfClass.containsKey(msn.getNameStringWithoutRt())){
//                  System.out.println("The positive list does not contain "+className+" "+msn.getNameStringWithoutRt()+"! But we identified it in "+fileName);
//                  continue;
//                }
//                LinkedHashMap<String,String> posStructures = posOfClass.get(msn.getNameStringWithoutRt());
//                boolean found = false;
//                if ((sheetName.equalsIgnoreCase("Orbitrap Velos Pro CID")||sheetName.equalsIgnoreCase("Orbitrap Velos Pro HCD"))
//                    && className.equalsIgnoreCase("DG")){
//                  faId = faId.replaceAll("/", "_");
//                  if (faId.split("_").length==2) faId += "_-";
//                }
//                for (String struct : posStructures.keySet()){
//                  if (StaticUtils.isAPermutedVersion(faId.replaceAll("/", "_"), struct.replaceAll("/", "_"))){
//                    found = true;
//                    break;
//                  }
//                }
//                if (!found) System.out.println("The positive list does not contain the structure "+className+msn.getNameStringWithoutRt()+" "+faId+"! But we identified it in "+fileName);
//              }else{
//                if (posOfClass==null || !posOfClass.containsKey(faId)) System.out.println("The positive list does not contain "+className+" "+msn.getNameStringWithoutRt()+"! But we identified it in "+fileName);
//              }
//            }
//          }
//        }
//      }
      
      
      
    }catch (Exception ex){
      ex.printStackTrace();
    }finally{
      if (out!=null){
        try {
          out.close();
        }
        catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
    }
    
  }
  
  private void generateDetailsExperiment1(){
    Vector<String> quantitationFiles = new Vector<String>();
    quantitationFiles.add("D:\\Experiment1\\massLists\\Ex1_pos.xlsx");
    quantitationFiles.add("D:\\Experiment1\\massLists\\Ex1_neg.xlsx");
    
    // these are the values for the G6550A QTOF
//    String resultFile = "D:\\Experiment1\\Singapore\\SpeciesDetectable_MSMS_G6550A_QTOF_Details_generated.xlsx";
//    String sheetName = "G6550A QTOF";
//    
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\positive\\Sample1_pos_01_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\positive\\Sample1_pos_02_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\positive\\Sample1_pos_03_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\positive\\Sample1_pos_04_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\positive\\Sample1_pos_05_Ex1_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//    
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\negative\\Sample1_neg_01_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\negative\\Sample1_neg_02_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\negative\\Sample1_neg_03_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\negative\\Sample1_neg_04_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Singapore\\negative\\Sample1_neg_05_Ex1_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.004d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getQTOFG6550Modifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getQTOFG6550PreElutionRanges();

    // these are the values for the Orbitrap Elite
//    String resultFile = "D:\\Experiment1\\Cambridge\\SpeciesDetectable_MSMS_Orbitrap_Elite_Details_generated.xlsx";
//    String sheetName = "Orbitrap Elite";
//    
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\positive\\2a_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\positive\\2b_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\positive\\2c_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\positive\\2d_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\positive\\2e_Ex1_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//    
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\negative\\2a_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\negative\\2b_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\negative\\2c_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\negative\\2d_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Cambridge\\negative\\2e_Ex1_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.003d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapEliteModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getOrbitrapElitePreElutionRanges();
    
    //these are the values for the Orbitrap_CID
//    String resultFile = "D:\\Experiment1\\Orbitrap_CID\\SpeciesDetectable_MSMS_Orbitrap_CID_Details_generated.xlsx";
//    String sheetName = "Orbitrap Velos Pro CID";
//    
//    
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\positive\\50\\014_Ex1_Orbitrap_CID_pos_50_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\positive\\50\\015_Ex1_Orbitrap_CID_pos_50_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\positive\\50\\016_Ex1_Orbitrap_CID_pos_50_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\positive\\50\\017_Ex1_Orbitrap_CID_pos_50_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\positive\\50\\018_Ex1_Orbitrap_CID_pos_50_Ex1_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//    
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\negative\\50\\014_Ex1_Orbitrap_CID_neg_50_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\negative\\50\\015_Ex1_Orbitrap_CID_neg_50_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\negative\\50\\016_Ex1_Orbitrap_CID_neg_50_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\negative\\50\\017_Ex1_Orbitrap_CID_neg_50_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_CID\\negative\\50\\018_Ex1_Orbitrap_CID_neg_50_Ex1_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//          
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.002d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapCIDModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getOrbitrapCIDPreElutionRanges();
    
    //these are the values for the Orbitrap_HCD
    String resultFile = "D:\\Experiment1\\Orbitrap_HCD\\SpeciesDetectable_MSMS_Orbitrap_HCD_Details_generated.xlsx";
    String sheetName = "Orbitrap Velos Pro HCD";
      
    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
    Vector<String> filesOfOneIonMode = new Vector<String>();
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\positive\\25\\044_Ex1_Orbitrap_HCD_pos_25_Ex1_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\positive\\25\\045_Ex1_Orbitrap_HCD_pos_25_Ex1_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\positive\\25\\046_Ex1_Orbitrap_HCD_pos_25_Ex1_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\positive\\25\\047_Ex1_Orbitrap_HCD_pos_25_Ex1_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\positive\\25\\048_Ex1_Orbitrap_HCD_pos_25_Ex1_pos.xlsx");
    ldaFilesIonMode.put("+", filesOfOneIonMode);
  
    filesOfOneIonMode = new Vector<String>();
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\negative\\50\\038_Ex1_Orbitrap_HCD_neg_50_Ex1_neg.xlsx");
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\negative\\50\\039_Ex1_Orbitrap_HCD_neg_50_Ex1_neg.xlsx");
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\negative\\50\\040_Ex1_Orbitrap_HCD_neg_50_Ex1_neg.xlsx");
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\negative\\50\\041_Ex1_Orbitrap_HCD_neg_50_Ex1_neg.xlsx");
    filesOfOneIonMode.add("D:\\Experiment1\\Orbitrap_HCD\\negative\\50\\042_Ex1_Orbitrap_HCD_neg_50_Ex1_neg.xlsx");
    ldaFilesIonMode.put("-", filesOfOneIonMode);
  
    int decimalPlacesForMz = 4;
    double highestMzStdevAllowed = 0.002d;
    double highestRtStdevAllowed = 0.2;
    double highestAbundanceCVAllowed = 50;
    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapHCDModifications();
    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getOrbitrapHCDPreElutionRanges();
    
    // these are the values for the Q Exactive
//    String resultFile = "D:\\Experiment1\\Stockholm\\SpeciesDetectable_MSMS_QExactive_Details_generated.xlsx";
//    String sheetName = "Q Exactive";
//      
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\positive\\20150709_LipidMix_5um-pos-01_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\positive\\20150709_LipidMix_5um-pos-02_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\positive\\20150709_LipidMix_5um-pos-03_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\positive\\20150709_LipidMix_5um-pos-04_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\positive\\20150709_LipidMix_5um-pos-05_Ex1_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//  
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\negative\\20150709_LipidMix_5um-neg-01_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\negative\\20150709_LipidMix_5um-neg-02_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\negative\\20150709_LipidMix_5um-neg-03_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\negative\\20150709_LipidMix_5um-neg-04_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\Stockholm\\negative\\20150709_LipidMix_5um-neg-05_Ex1_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//  
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.006d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getQExactiveModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getQExactivePreElutionRanges();

    // these are the values for the 4000 QTRAP
//    String resultFile = "D:\\Experiment1\\QTRAP\\SpeciesDetectable_MSMS_4000QTRAP_Details_generated.xlsx";
//    String sheetName = "4000 QTRAP";
//      
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\positive\\45\\Data20150309_Ex1_pos-014_Ex1_Qtrap_pos_45_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\positive\\45\\Data20150309_Ex1_pos-015_Ex1_Qtrap_pos_45_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\positive\\45\\Data20150309_Ex1_pos-016_Ex1_Qtrap_pos_45_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\positive\\45\\Data20150309_Ex1_pos-017_Ex1_Qtrap_pos_45_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\positive\\45\\Data20150309_Ex1_pos-018_Ex1_Qtrap_pos_45_Ex1_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//  
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\negative\\45\\Data20150105_Ex1_neg-014_Ex1_Qtrap_neg_45_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\negative\\45\\Data20150105_Ex1_neg-015_Ex1_Qtrap_neg_45_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\negative\\45\\Data20150105_Ex1_neg-016_Ex1_Qtrap_neg_45_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\negative\\45\\Data20150105_Ex1_neg-017_Ex1_Qtrap_neg_45_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTRAP\\negative\\45\\Data20150105_Ex1_neg-018_Ex1_Qtrap_neg_45_Ex1_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//  
//    int decimalPlacesForMz = 2;
//    double highestMzStdevAllowed = 0.15d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.get4000QTRAPModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.get4000QTRAPPreElutionRanges();

    // these are the values for the QTRAP 6500
//    String resultFile = "D:\\Experiment1\\SanDiego\\SpeciesDetectable_MSMS_QTRAP6500_Details_generated.xlsx";
//    String sheetName = "QTRAP 6500";
//      
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\SanDiego\\positive\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Standard 1 5uM pos_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\SanDiego\\positive\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Standard 1 5uM pos (2)_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\SanDiego\\positive\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Standard 1 5uM pos (3)_Ex1_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//  
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\SanDiego\\negative\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Standard 1 5uM neg_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\SanDiego\\negative\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Standard 1 5uM neg (2)_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\SanDiego\\negative\\10022015 Hartler (4 Samples) EMS IDA neg and pos mode 5x replicates-Standard 1 5uM neg (3)_Ex1_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//  
//    int decimalPlacesForMz = 2;
//    double highestMzStdevAllowed = 0.15d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getQTRAP6500Modifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getQTRAP6500PreElutionRanges();

    
    // these are the values for the SYNAPT G1 HDMS QTOF
//    String resultFile = "D:\\Experiment1\\QTOF\\SpeciesDetectable_MSMS_SynaptG1_Details_generated.xlsx";
//    String sheetName = "SYNAPT G1 HDMS QTOF";
//      
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\positive\\30\\20150825_GNR_LDApos_Coll30_01_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\positive\\30\\20150825_GNR_LDApos_Coll30_02_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\positive\\30\\20150825_GNR_LDApos_Coll30_03_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\positive\\30\\20150825_GNR_LDApos_Coll30_04_Ex1_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\positive\\30\\20150825_GNR_LDApos_Coll30_05_Ex1_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//  
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\negative\\30\\201500706_GNR_FragTest_CollEng_30neg_01_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\negative\\30\\201500706_GNR_FragTest_CollEng_30neg_02_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\negative\\30\\201500706_GNR_FragTest_CollEng_30neg_03_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\negative\\30\\201500706_GNR_FragTest_CollEng_30neg_04_Ex1_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment1\\QTOF\\negative\\30\\201500706_GNR_FragTest_CollEng_30neg_05_Ex1_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//  
//    int decimalPlacesForMz = 3;
//    double highestMzStdevAllowed = 0.015d;
//    double highestRtStdevAllowed = 1.0;
//    double highestAbundanceCVAllowed = 75;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getSynaptG1Modifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getSynaptG1PreElutionRanges();

    
    BufferedOutputStream out = null;
    Hashtable<String,String> appliedMzCorrection = new Hashtable<String,String>();
    int highestNumberOfFiles = 0;
    if (ldaFilesIonMode.get("+").size()>highestNumberOfFiles) highestNumberOfFiles = ldaFilesIonMode.get("+").size();
    if (ldaFilesIonMode.get("-").size()>highestNumberOfFiles) highestNumberOfFiles = ldaFilesIonMode.get("-").size();

    try{
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> theoreticalMasses = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();          
      for (String quantFile : quantitationFiles){
        Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> oneFile = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)QuantificationThread.parseQuantExcelFile(quantFile,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f,false).get(3);
        ////Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> oneFile = null;
        for (String className : oneFile.keySet()){
          Hashtable<String,Hashtable<String,QuantVO>> oneClass= oneFile.get(className);
          Hashtable<String,Hashtable<String,QuantVO>> ofClass = new Hashtable<String,Hashtable<String,QuantVO>>();
          if (theoreticalMasses.containsKey(className)) ofClass = theoreticalMasses.get(className);
          for (String analyte:oneClass.keySet()){
            Hashtable<String,QuantVO> oneAnalyte = oneClass.get(analyte);
            String analyteString = new String(analyte);
            if (className.equalsIgnoreCase("TG")&&analyteString.startsWith("d")) analyteString = analyteString.substring(1);
            Hashtable<String,QuantVO> ofAnalyte = new Hashtable<String,QuantVO>();
            if (ofClass.containsKey(analyteString)) ofAnalyte = ofClass.get(analyteString);
            for (String mod:oneAnalyte.keySet()) ofAnalyte.put(mod, oneAnalyte.get(mod));
            ofClass.put(analyteString, ofAnalyte);
          }
          theoreticalMasses.put(className,ofClass);
        }
      }
   
    
    
      LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> positiveList = this.getPositiveListOfExp1Standards();
//      for (String className : positiveList.keySet()){
//        System.out.println(className);
//        for (String species : positiveList.get(className).keySet()){
//          System.out.println("     "+species);
//          for (String structure : positiveList.get(className).get(species).values()){
//            System.out.println("            "+structure);
//          }
//        }
//      }
      
      //this hashtable is for checking if all of the species are in the positive list
      LinkedHashMap<String,LinkedHashMap<String,Vector<LipidParameterSet>>> checkPositiveList = new LinkedHashMap<String,LinkedHashMap<String,Vector<LipidParameterSet>>>();
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>> results = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>>();
      for (String polarity : ldaFilesIonMode.keySet()){
        for (String ldaFile : ldaFilesIonMode.get(polarity)){
          String fileName = ldaFile.substring(ldaFile.lastIndexOf("\\")+1);
          QuantificationResult quantResult = LDAResultReader.readResultFile(ldaFile,  new Hashtable<String,Boolean>());
          if (quantResult.getConstants().getShift()==0d) appliedMzCorrection.put(fileName, "0");
          else{
            System.out.println(quantResult.getConstants().getShift());
            appliedMzCorrection.put(fileName, Calculator.FormatNumberToString(quantResult.getConstants().getShift(),decimalPlacesForMz));
          }
          Hashtable<String,Vector<LipidParameterSet>> idents = quantResult.getIdentifications();
          //this is for checking if all classes are present
          for (String className : idents.keySet()){
            Vector<LipidParameterSet> classResult = idents.get(className);
            //this is for checking if all classes are present
            boolean found = false;
            for (LipidParameterSet set : classResult){
              if (set instanceof LipidomicsMSnSet){
                found = true;
                break;
              }
            }
            if (found && !positiveList.containsKey(className)){
              System.out.println("The positive list does not contain the lipid class "+className+"! But we identified it in "+fileName);
            }
          // end of checking if all classes are present
          
          //this is for checking if all of the species are in the positive list
//            LinkedHashMap<String,Vector<LipidParameterSet>> classCheck = new LinkedHashMap<String,Vector<LipidParameterSet>>();
//            if (checkPositiveList.containsKey(className)) classCheck = checkPositiveList.get(className);
//            classCheck.put(fileName, classResult);
//            checkPositiveList.put(className, classCheck);
            // end of checking if all of the species are in the positive list

          }
        
          Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> resultsOfExp = new Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>();
          for (String className : positiveList.keySet()){
            if (!idents.containsKey(className)) continue;
            Vector<LipidParameterSet> classResult = idents.get(className);
            LinkedHashMap<String,LinkedHashMap<String,String>> posAnalytes = positiveList.get(className);
            Hashtable<String,Vector<LipidomicsMSnSet>> msnFound = new Hashtable<String,Vector<LipidomicsMSnSet>>();
            for (LipidParameterSet set : classResult){
              String analyte = set.getNameStringWithoutRt();
              if (className.equalsIgnoreCase("TG") && analyte.startsWith("d")) analyte = analyte.substring(1);
              if (!(set instanceof LipidomicsMSnSet) || !posAnalytes.containsKey(analyte)) continue;
              float rt = Float.parseFloat(set.getRt());
              if (sheetName.equalsIgnoreCase("Q Exactive") && (className.equalsIgnoreCase("TG"))){
                if (set.getNameStringWithoutRt().equalsIgnoreCase("58:5") && rt<41.7f) continue;
                else if (set.getNameStringWithoutRt().equalsIgnoreCase("58:7") && rt<40.45f) continue;
                
              }
                  
              Vector<LipidomicsMSnSet> ofOneAnalyte = new Vector<LipidomicsMSnSet>();
              if (msnFound.containsKey(analyte)) ofOneAnalyte = msnFound.get(analyte);
                ofOneAnalyte.add((LipidomicsMSnSet)set);
                msnFound.put(analyte, ofOneAnalyte);
            }
            resultsOfExp.put(className, msnFound);
          }
          results.put(fileName, resultsOfExp);
        }
      }
      
      Workbook resultWorkbook = createDetailsExcelFile(sheetName,highestNumberOfFiles,decimalPlacesForMz,highestMzStdevAllowed,highestRtStdevAllowed,
          highestAbundanceCVAllowed,ldaFilesIonMode,appliedMzCorrection,mods,preElutionRanges,positiveList,theoreticalMasses,results,true,false);
      out = new BufferedOutputStream(new FileOutputStream(resultFile));
      resultWorkbook.write(out);
      resultWorkbook.close();

      
      // this is for checking the positive list
//    for (String className : checkPositiveList.keySet()){
//      LinkedHashMap<String,Vector<LipidParameterSet>> result = checkPositiveList.get(className);
//      for (String fileName : result.keySet()){
//        Vector<LipidParameterSet> classResults = result.get(fileName);
//        for (LipidParameterSet set : classResults){
//          if (!(set instanceof LipidomicsMSnSet)) continue;
//          LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
//          for (Object msnNames : msn.getMSnIdentificationNames()){    
//            String nameString = "";
//            String faId = "";
//            if (msnNames instanceof Vector){
//              for (String name : (Vector<String>)msnNames){
//                nameString += name+"|";
//              }
//              faId = nameString.substring(0,nameString.indexOf("|"));
//              nameString = nameString.substring(0,nameString.length()-1);
//            }else{
//              nameString = (String)msnNames;
//              faId = nameString;
//            }
//            LinkedHashMap<String,LinkedHashMap<String,String>> posOfClass = positiveList.get(className);
//            if (msn.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED && !faId.equalsIgnoreCase(msn.getNameStringWithoutRt())){
//              String analyte = set.getNameStringWithoutRt();
//              if (className.equalsIgnoreCase("TG") && analyte.startsWith("d")) analyte = analyte.substring(1);
//              if (!posOfClass.containsKey(analyte)){
//                System.out.println("The positive list does not contain "+className+" "+analyte+"! But we identified it in "+fileName);
//                continue;
//              }
//              LinkedHashMap<String,String> posStructures = posOfClass.get(analyte);
//              boolean found = false;
//              if (//(sheetName.equalsIgnoreCase("Orbitrap Velos Pro CID")||sheetName.equalsIgnoreCase("Orbitrap Velos Pro HCD"))&& 
//                  className.equalsIgnoreCase("DG")){
//                faId = faId.replaceAll("/", "_");
//                if (faId.split("_").length==2) faId += "_-";
//              }
//              for (String struct : posStructures.keySet()){
//                if (StaticUtils.isAPermutedVersion(faId.replaceAll("/", "_"), struct.replaceAll("/", "_"))){
//                  found = true;
//                  break;
//                }
//              }
//              if (!found) System.out.println("The positive list does not contain the structure "+className+analyte+" "+faId+"! But we identified it in "+fileName);
//            }else{
//              String analyte = set.getNameStringWithoutRt();
//              if (className.equalsIgnoreCase("TG") && analyte.startsWith("d")) analyte = analyte.substring(1);
//              if (posOfClass==null || (!posOfClass.containsKey(analyte))) System.out.println("The positive list does not contain "+className+" "+msn.getNameStringWithoutRt()+"! But we identified it in "+fileName);
//            }
//          }
//        }
//      }
//    }
      
      
    }catch (Exception ex){
      ex.printStackTrace();
    }finally{
      if (out!=null){
        try {
          out.close();
        }
        catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
    }
  }

  private void generateDetailsExperiment2(){
    Vector<String> quantitationFiles = new Vector<String>();
    quantitationFiles.add("D:\\Experiment2\\massLists\\Ex2_pos.xlsx");
    
    
    //these are the values for the Orbitrap_CID
    String resultFile = "D:\\Experiment2\\Orbitrap_CID\\Experiment2_MSMS_Orbitrap_CID_Details_generated.xlsx";
    String sheetName = "Orbitrap Velos Pro CID";
  
    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
    Vector<String> filesOfOneIonMode = new Vector<String>();
    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_CID\\001_Ex2_Orbitrap_CID_Ex2_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_CID\\002_Ex2_Orbitrap_CID_Ex2_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_CID\\003_Ex2_Orbitrap_CID_Ex2_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_CID\\004_Ex2_Orbitrap_CID_Ex2_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_CID\\005_Ex2_Orbitrap_CID_Ex2_pos.xlsx");
    ldaFilesIonMode.put("+", filesOfOneIonMode);
  
    int decimalPlacesForMz = 4;
    double highestMzStdevAllowed = 0.002d;
    double highestRtStdevAllowed = 0.2;
    double highestAbundanceCVAllowed = 50;
    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapCIDModifications();
    Hashtable<String,Range> preElutionRanges = new Hashtable<String,Range>();

    
    //these are the values for the Orbitrap_HCD
//    String resultFile = "D:\\Experiment2\\Orbitrap_HCD\\Experiment2_MSMS_Orbitrap_HCD_Details_generated.xlsx";
//    String sheetName = "Orbitrap Velos Pro HCD";
//    
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_HCD\\001_Ex2_Orbitrap_HCD_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_HCD\\002_Ex2_Orbitrap_HCD_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_HCD\\003_Ex2_Orbitrap_HCD_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_HCD\\004_Ex2_Orbitrap_HCD_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\Orbitrap_HCD\\005_Ex2_Orbitrap_HCD_Ex2_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.002d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapHCDModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getOrbitrapHCDPreElutionRanges();

    // these are the values for the 4000 QTRAP
//    String resultFile = "D:\\Experiment2\\QTRAP\\Experiment2_MSMS_4000QTRAP_Details_generated.xlsx";
//    String sheetName = "4000 QTRAP";
//    
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment2\\QTRAP\\Ex2_QTrap93-001_Ex2_QTrap_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\QTRAP\\Ex2_QTrap93-002_Ex2_QTrap_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\QTRAP\\Ex2_QTrap93-003_Ex2_QTrap_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\QTRAP\\Ex2_QTrap93-004_Ex2_QTrap_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\QTRAP\\Ex2_QTrap93-005_Ex2_QTrap_Ex2_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//
//    int decimalPlacesForMz = 2;
//    double highestMzStdevAllowed = 0.15d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.get4000QTRAPModifications();
 //   Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.get4000QTRAPPreElutionRanges();

    // these are the values for the SYNAPT G1 HDMS QTOF
//    String resultFile = "D:\\Experiment2\\QTOF\\Experiment2_MSMS_SynaptG1_Details_generated.xlsx";
//    String sheetName = "SYNAPT G1 HDMS QTOF";
//      
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment2\\QTOF\\20151028_GNR_LDA_Exp2_01_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\QTOF\\20151028_GNR_LDA_Exp2_02_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\QTOF\\20151028_GNR_LDA_Exp2_03_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\QTOF\\20151028_GNR_LDA_Exp2_04_Ex2_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment2\\QTOF\\20151028_GNR_LDA_Exp2_05_Ex2_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//    
//    int decimalPlacesForMz = 3;
//    double highestMzStdevAllowed = 0.015d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 55;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getSynaptG1Modifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getSynaptG1PreElutionRanges();

    BufferedOutputStream out = null;
    Hashtable<String,String> appliedMzCorrection = new Hashtable<String,String>();
    int highestNumberOfFiles = 0;
    if (ldaFilesIonMode.get("+").size()>highestNumberOfFiles) highestNumberOfFiles = ldaFilesIonMode.get("+").size();

    try{
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> theoreticalMasses = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();          
      for (String quantFile : quantitationFiles){
        Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> oneFile = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)QuantificationThread.parseQuantExcelFile(quantFile,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f,false).get(3);
        ////Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> oneFile = null;
        for (String className : oneFile.keySet()){
          Hashtable<String,Hashtable<String,QuantVO>> oneClass= oneFile.get(className);
          Hashtable<String,Hashtable<String,QuantVO>> ofClass = new Hashtable<String,Hashtable<String,QuantVO>>();
          if (theoreticalMasses.containsKey(className)) ofClass = theoreticalMasses.get(className);
          for (String analyte:oneClass.keySet()){
            Hashtable<String,QuantVO> oneAnalyte = oneClass.get(analyte);
            String analyteString = new String(analyte);
            Hashtable<String,QuantVO> ofAnalyte = new Hashtable<String,QuantVO>();
            if (ofClass.containsKey(analyteString)) ofAnalyte = ofClass.get(analyteString);
            for (String mod:oneAnalyte.keySet()) ofAnalyte.put(mod, oneAnalyte.get(mod));
            ofClass.put(analyteString, ofAnalyte);
          }
          theoreticalMasses.put(className,ofClass);
        }
      }
      
      LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> positiveList = this.getPositiveListOfExp2Standards();
//    for (String className : positiveList.keySet()){
//      System.out.println(className);
//      for (String species : positiveList.get(className).keySet()){
//        System.out.println("     "+species);
//        for (String structure : positiveList.get(className).get(species).values()){
//          System.out.println("            "+structure);
//        }
//      }
//    }
      //this hashtable is for checking if all of the species are in the positive list
      LinkedHashMap<String,LinkedHashMap<String,Vector<LipidParameterSet>>> checkPositiveList = new LinkedHashMap<String,LinkedHashMap<String,Vector<LipidParameterSet>>>();
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>> results = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>>();
      for (String polarity : ldaFilesIonMode.keySet()){
        for (String ldaFile : ldaFilesIonMode.get(polarity)){
          String fileName = ldaFile.substring(ldaFile.lastIndexOf("\\")+1);
          QuantificationResult quantResult = LDAResultReader.readResultFile(ldaFile,  new Hashtable<String,Boolean>());
          if (quantResult.getConstants().getShift()==0d) appliedMzCorrection.put(fileName, "0");
          else{
            System.out.println(quantResult.getConstants().getShift());
            appliedMzCorrection.put(fileName, Calculator.FormatNumberToString(quantResult.getConstants().getShift(),decimalPlacesForMz));
          }
          Hashtable<String,Vector<LipidParameterSet>> idents = quantResult.getIdentifications();
          //this is for checking if all classes are present
          for (String className : idents.keySet()){
            Vector<LipidParameterSet> classResult = idents.get(className);
            //this is for checking if all classes are present
            boolean found = false;
            for (LipidParameterSet set : classResult){
              if (set instanceof LipidomicsMSnSet){
                found = true;
                break;
              }
            }
            if (found && !positiveList.containsKey(className)){
              System.out.println("The positive list does not contain the lipid class "+className+"! But we identified it in "+fileName);
            }
          // end of checking if all classes are present
          
          //this is for checking if all of the species are in the positive list
//            LinkedHashMap<String,Vector<LipidParameterSet>> classCheck = new LinkedHashMap<String,Vector<LipidParameterSet>>();
//            if (checkPositiveList.containsKey(className)) classCheck = checkPositiveList.get(className);
//            classCheck.put(fileName, classResult);
//            checkPositiveList.put(className, classCheck);
            // end of checking if all of the species are in the positive list

          }
        
          Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> resultsOfExp = new Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>();
          for (String className : positiveList.keySet()){
            if (!idents.containsKey(className)) continue;
            Vector<LipidParameterSet> classResult = idents.get(className);
            LinkedHashMap<String,LinkedHashMap<String,String>> posAnalytes = positiveList.get(className);
            Hashtable<String,Vector<LipidomicsMSnSet>> msnFound = new Hashtable<String,Vector<LipidomicsMSnSet>>();
            for (LipidParameterSet set : classResult){
              String analyte = set.getNameStringWithoutRt();
              if (!(set instanceof LipidomicsMSnSet) || !posAnalytes.containsKey(analyte)) continue;
              float rt = Float.parseFloat(set.getRt());
              Vector<LipidomicsMSnSet> ofOneAnalyte = new Vector<LipidomicsMSnSet>();
              if (msnFound.containsKey(analyte)) ofOneAnalyte = msnFound.get(analyte);
                ofOneAnalyte.add((LipidomicsMSnSet)set);
                msnFound.put(analyte, ofOneAnalyte);
            }
            resultsOfExp.put(className, msnFound);
          }
          results.put(fileName, resultsOfExp);
        }
      }

      Workbook resultWorkbook = createDetailsExcelFile(sheetName,highestNumberOfFiles,decimalPlacesForMz,highestMzStdevAllowed,highestRtStdevAllowed,
          highestAbundanceCVAllowed,ldaFilesIonMode,appliedMzCorrection,mods,preElutionRanges,positiveList,theoreticalMasses,results,true,false);
      out = new BufferedOutputStream(new FileOutputStream(resultFile));
      resultWorkbook.write(out);
      resultWorkbook.close();

      
      // this is for checking the positive list
//    for (String className : checkPositiveList.keySet()){
//      LinkedHashMap<String,Vector<LipidParameterSet>> result = checkPositiveList.get(className);
//      for (String fileName : result.keySet()){
//        Vector<LipidParameterSet> classResults = result.get(fileName);
//        for (LipidParameterSet set : classResults){
//          if (!(set instanceof LipidomicsMSnSet)) continue;
//          LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
//          for (Object msnNames : msn.getMSnIdentificationNames()){    
//            String nameString = "";
//            String faId = "";
//            if (msnNames instanceof Vector){
//              for (String name : (Vector<String>)msnNames){
//                nameString += name+"|";
//              }
//              faId = nameString.substring(0,nameString.indexOf("|"));
//              nameString = nameString.substring(0,nameString.length()-1);
//            }else{
//              nameString = (String)msnNames;
//              faId = nameString;
//            }
//            LinkedHashMap<String,LinkedHashMap<String,String>> posOfClass = positiveList.get(className);
//            if (msn.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED && !faId.equalsIgnoreCase(msn.getNameStringWithoutRt())){
//              String analyte = set.getNameStringWithoutRt();
//              if (className.equalsIgnoreCase("TG") && analyte.startsWith("d")) analyte = analyte.substring(1);
//              if (!posOfClass.containsKey(analyte)){
//                System.out.println("The positive list does not contain "+className+" "+analyte+"! But we identified it in "+fileName);
//                continue;
//              }
//              LinkedHashMap<String,String> posStructures = posOfClass.get(analyte);
//              boolean found = false;
//              if (//(sheetName.equalsIgnoreCase("Orbitrap Velos Pro CID")||sheetName.equalsIgnoreCase("Orbitrap Velos Pro HCD"))&& 
//                  className.equalsIgnoreCase("DG")){
//                faId = faId.replaceAll("/", "_");
//                if (faId.split("_").length==2) faId += "_-";
//              }
//              for (String struct : posStructures.keySet()){
//                if (StaticUtils.isAPermutedVersion(faId.replaceAll("/", "_"), struct.replaceAll("/", "_"))){
//                  found = true;
//                  break;
//                }
//              }
//              if (!found) System.out.println("The positive list does not contain the structure "+className+analyte+" "+faId+"! But we identified it in "+fileName);
//            }else{
//              String analyte = set.getNameStringWithoutRt();
//              if (className.equalsIgnoreCase("TG") && analyte.startsWith("d")) analyte = analyte.substring(1);
//              if (posOfClass==null || (!posOfClass.containsKey(analyte))) System.out.println("The positive list does not contain "+className+" "+msn.getNameStringWithoutRt()+"! But we identified it in "+fileName);
//            }
//          }
//        }
//      }
//    }
      
      
    }catch (Exception ex){
      ex.printStackTrace();
    }finally{
      if (out!=null){
        try {
          out.close();
        }
        catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
    }
  }
  
  private void generateDetailsExperiment3(){
    Vector<String> quantitationFiles = new Vector<String>();
    quantitationFiles.add("D:\\Experiment3\\massLists\\Ex3_pos.xlsx");
    quantitationFiles.add("D:\\Experiment3\\massLists\\Ex3_neg.xlsx");

    //these are the values for the Orbitrap_CID
//    String resultFile = "D:\\Experiment3\\Orbitrap_CID\\Experiment3_MSMS_Orbitrap_CID_Details_generated.xlsx";
//    String sheetName = "Orbitrap Velos Pro CID";
//  
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\positive\\002_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\positive\\003_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\positive\\004_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\positive\\005_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\positive\\006_Ex3-1_Orbitrap_CID_pos_Ex3_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//  
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\negative\\002_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\negative\\003_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\negative\\004_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\negative\\005_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_CID\\negative\\006_Ex3-1_Orbitrap_CID_neg_Ex3_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//        
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.002d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapCIDModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getOrbitrapCIDPreElutionRanges();

    
    // these are the values for the Orbitrap_HCD
//    String resultFile = "D:\\Experiment3\\Orbitrap_HCD\\Experiment3_MSMS_Orbitrap_HCD_Details_generated.xlsx";
//    String sheetName = "Orbitrap Velos Pro HCD";
//      
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\positive\\034_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\positive\\035_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\positive\\036_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\positive\\037_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\positive\\038_Ex3-1_Orbitrap_HCD_pos_Ex3_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//  
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\negative\\034_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\negative\\035_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\negative\\036_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\negative\\037_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\Orbitrap_HCD\\negative\\038_Ex3-1_Orbitrap_HCD_neg_Ex3_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//  
//    int decimalPlacesForMz = 4;
//    double highestMzStdevAllowed = 0.002d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getOrbitrapHCDModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getOrbitrapHCDPreElutionRanges();

    // these are the values for the 4000 QTRAP
//    String resultFile = "D:\\Experiment3\\QTRAP\\Experiment3_MSMS_4000QTRAP_Details_generated.xlsx";
//    String sheetName = "4000 QTRAP";
//    
//    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
//    Vector<String> filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\positive\\Data20150827_Ex3_QTrap_pos-001_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\positive\\Data20150827_Ex3_QTrap_pos-002_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\positive\\Data20150827_Ex3_QTrap_pos-003_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\positive\\Data20150827_Ex3_QTrap_pos-004_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\positive\\Data20150827_Ex3_QTrap_pos-005_QTrap_Ex3.1_pos_Ex3_pos.xlsx");
//    ldaFilesIonMode.put("+", filesOfOneIonMode);
//
//    filesOfOneIonMode = new Vector<String>();
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\negative\\Data20150826_Ex3_QTrap_neg-001_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\negative\\Data20150826_Ex3_QTrap_neg-002_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\negative\\Data20150826_Ex3_QTrap_neg-003_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\negative\\Data20150826_Ex3_QTrap_neg-004_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//    filesOfOneIonMode.add("D:\\Experiment3\\QTRAP\\negative\\Data20150826_Ex3_QTrap_neg-005_QTrap_Ex3.1_neg_Ex3_neg.xlsx");
//    ldaFilesIonMode.put("-", filesOfOneIonMode);
//
//    int decimalPlacesForMz = 2;
//    double highestMzStdevAllowed = 0.15d;
//    double highestRtStdevAllowed = 0.2;
//    double highestAbundanceCVAllowed = 50;
//    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.get4000QTRAPModifications();
//    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.get4000QTRAPPreElutionRanges();

    // these are the values for the SYNAPT G1 HDMS QTOF
    String resultFile = "D:\\Experiment3\\QTOF\\Experiment3_MSMS_SynaptG1_Details_generated.xlsx";
    String sheetName = "SYNAPT G1 HDMS QTOF";
    
    LinkedHashMap<String,Vector<String>> ldaFilesIonMode = new LinkedHashMap<String,Vector<String>>();
    Vector<String> filesOfOneIonMode = new Vector<String>();
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\positive\\20151103_GNR_LDA_Exp3_1_01_Ex3_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\positive\\20151103_GNR_LDA_Exp3_1_02_Ex3_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\positive\\20151103_GNR_LDA_Exp3_1_03_Ex3_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\positive\\20151103_GNR_LDA_Exp3_1_04_Ex3_pos.xlsx");
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\positive\\20151103_GNR_LDA_Exp3_1_05_Ex3_pos.xlsx");
    ldaFilesIonMode.put("+", filesOfOneIonMode);

    filesOfOneIonMode = new Vector<String>();
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\negative\\20151103_GNR_LDA_Exp3_1_01n_Ex3_neg.xlsx");
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\negative\\20151103_GNR_LDA_Exp3_1_02n_Ex3_neg.xlsx");
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\negative\\20151103_GNR_LDA_Exp3_1_03n_Ex3_neg.xlsx");
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\negative\\20151103_GNR_LDA_Exp3_1_04n_Ex3_neg.xlsx");
    filesOfOneIonMode.add("D:\\Experiment3\\QTOF\\negative\\20151103_GNR_LDA_Exp3_1_05n_Ex3_neg.xlsx");
    ldaFilesIonMode.put("-", filesOfOneIonMode);

    int decimalPlacesForMz = 3;
    double highestMzStdevAllowed = 0.015d;
    double highestRtStdevAllowed = 0.2;
    double highestAbundanceCVAllowed = 55;
    Hashtable<String,LinkedHashMap<String,String>> mods = FoundBiologicalSpecies.getSynaptG1Modifications();
    Hashtable<String,Range> preElutionRanges = FoundBiologicalSpecies.getSynaptG1PreElutionRanges();

    
    BufferedOutputStream out = null;
    Hashtable<String,String> appliedMzCorrection = new Hashtable<String,String>();
    int highestNumberOfFiles = 0;
    if (ldaFilesIonMode.get("+").size()>highestNumberOfFiles) highestNumberOfFiles = ldaFilesIonMode.get("+").size();
    if (ldaFilesIonMode.get("-").size()>highestNumberOfFiles) highestNumberOfFiles = ldaFilesIonMode.get("-").size();

    try{
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> theoreticalMasses = new Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>();          
      for (String quantFile : quantitationFiles){
        Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> oneFile = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)QuantificationThread.parseQuantExcelFile(quantFile,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f,false).get(3);
        ////Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> oneFile = null;
        for (String className : oneFile.keySet()){
          Hashtable<String,Hashtable<String,QuantVO>> oneClass= oneFile.get(className);
          Hashtable<String,Hashtable<String,QuantVO>> ofClass = new Hashtable<String,Hashtable<String,QuantVO>>();
          if (theoreticalMasses.containsKey(className)) ofClass = theoreticalMasses.get(className);
          for (String analyte:oneClass.keySet()){
            Hashtable<String,QuantVO> oneAnalyte = oneClass.get(analyte);
            String analyteString = new String(analyte);
            Hashtable<String,QuantVO> ofAnalyte = new Hashtable<String,QuantVO>();
            if (ofClass.containsKey(analyteString)) ofAnalyte = ofClass.get(analyteString);
            for (String mod:oneAnalyte.keySet()) ofAnalyte.put(mod, oneAnalyte.get(mod));
            ofClass.put(analyteString, ofAnalyte);
          }
          theoreticalMasses.put(className,ofClass);
        }
      }
      LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> positiveList = this.getPositiveListOfExp3Standards();
//      for (String className : positiveList.keySet()){
//        System.out.println(className);
//        for (String species : positiveList.get(className).keySet()){
//          System.out.println("     "+species);
//          for (String structure : positiveList.get(className).get(species).values()){
//            System.out.println("            "+structure);
//          }
//        }
//      }

      //this hashtable is for checking if all of the species are in the positive list
      LinkedHashMap<String,LinkedHashMap<String,Vector<LipidParameterSet>>> checkPositiveList = new LinkedHashMap<String,LinkedHashMap<String,Vector<LipidParameterSet>>>();
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>> results = new Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>>();
      for (String polarity : ldaFilesIonMode.keySet()){
        for (String ldaFile : ldaFilesIonMode.get(polarity)){
          String fileName = ldaFile.substring(ldaFile.lastIndexOf("\\")+1);
          QuantificationResult quantResult = LDAResultReader.readResultFile(ldaFile,  new Hashtable<String,Boolean>());
          if (quantResult.getConstants().getShift()==0d) appliedMzCorrection.put(fileName, "0");
          else{
            System.out.println(quantResult.getConstants().getShift());
            appliedMzCorrection.put(fileName, Calculator.FormatNumberToString(quantResult.getConstants().getShift(),decimalPlacesForMz));
          }
          Hashtable<String,Vector<LipidParameterSet>> idents = quantResult.getIdentifications();
          //this is for checking if all classes are present
          for (String className : idents.keySet()){
            Vector<LipidParameterSet> classResult = idents.get(className);
            //this is for checking if all classes are present
            boolean found = false;
            for (LipidParameterSet set : classResult){
              if (set instanceof LipidomicsMSnSet){
                found = true;
////                break;
              }
            }
            if (found && !positiveList.containsKey(className)){
              System.out.println("The positive list does not contain the lipid class "+className+"! But we identified it in "+fileName);
            }
          // end of checking if all classes are present
          
          //this is for checking if all of the species are in the positive list
//            LinkedHashMap<String,Vector<LipidParameterSet>> classCheck = new LinkedHashMap<String,Vector<LipidParameterSet>>();
//            if (checkPositiveList.containsKey(className)) classCheck = checkPositiveList.get(className);
//            classCheck.put(fileName, classResult);
//            checkPositiveList.put(className, classCheck);
            // end of checking if all of the species are in the positive list

          }
        
          Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> resultsOfExp = new Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>();
          for (String className : positiveList.keySet()){
            if (!idents.containsKey(className)) continue;
            Vector<LipidParameterSet> classResult = idents.get(className);
            LinkedHashMap<String,LinkedHashMap<String,String>> posAnalytes = positiveList.get(className);
            Hashtable<String,Vector<LipidomicsMSnSet>> msnFound = new Hashtable<String,Vector<LipidomicsMSnSet>>();
            for (LipidParameterSet set : classResult){
              String analyte = set.getNameStringWithoutRt();
              if (!(set instanceof LipidomicsMSnSet) || !posAnalytes.containsKey(analyte)) continue;
              float rt = Float.parseFloat(set.getRt());
              Vector<LipidomicsMSnSet> ofOneAnalyte = new Vector<LipidomicsMSnSet>();
              if (msnFound.containsKey(analyte)) ofOneAnalyte = msnFound.get(analyte);
                ofOneAnalyte.add((LipidomicsMSnSet)set);
                msnFound.put(analyte, ofOneAnalyte);
            }
            resultsOfExp.put(className, msnFound);
          }
          results.put(fileName, resultsOfExp);
        }
      }

      Workbook resultWorkbook = createDetailsExcelFile(sheetName,highestNumberOfFiles,decimalPlacesForMz,highestMzStdevAllowed,highestRtStdevAllowed,
          highestAbundanceCVAllowed,ldaFilesIonMode,appliedMzCorrection,mods,preElutionRanges,positiveList,theoreticalMasses,results,false,true);
      out = new BufferedOutputStream(new FileOutputStream(resultFile));
      resultWorkbook.write(out);
      resultWorkbook.close();

      
      // this is for checking the positive list
//    for (String className : checkPositiveList.keySet()){
//      LinkedHashMap<String,Vector<LipidParameterSet>> result = checkPositiveList.get(className);
//      for (String fileName : result.keySet()){
//        Vector<LipidParameterSet> classResults = result.get(fileName);
//        for (LipidParameterSet set : classResults){
//          if (!(set instanceof LipidomicsMSnSet)) continue;
//          LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
//          for (Object msnNames : msn.getMSnIdentificationNames()){    
//            String nameString = "";
//            String faId = "";
//            if (msnNames instanceof Vector){
//              for (String name : (Vector<String>)msnNames){
//                nameString += name+"|";
//              }
//              faId = nameString.substring(0,nameString.indexOf("|"));
//              nameString = nameString.substring(0,nameString.length()-1);
//            }else{
//              nameString = (String)msnNames;
//              faId = nameString;
//            }
//            LinkedHashMap<String,LinkedHashMap<String,String>> posOfClass = positiveList.get(className);
//            if (msn.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED && !faId.equalsIgnoreCase(msn.getNameStringWithoutRt())){
//              String analyte = set.getNameStringWithoutRt();
//              if (className.equalsIgnoreCase("TG") && analyte.startsWith("d")) analyte = analyte.substring(1);
//              if (!posOfClass.containsKey(analyte)){
//                System.out.println("The positive list does not contain "+className+" "+analyte+"! But we identified it in "+fileName);
//                continue;
//              }
//              LinkedHashMap<String,String> posStructures = posOfClass.get(analyte);
//              boolean found = false;
//              if (//(sheetName.equalsIgnoreCase("Orbitrap Velos Pro CID")||sheetName.equalsIgnoreCase("Orbitrap Velos Pro HCD"))&& 
//                  className.equalsIgnoreCase("DG")){
//                faId = faId.replaceAll("/", "_");
//                if (faId.split("_").length==2) faId += "_-";
//              }
//              for (String struct : posStructures.keySet()){
//                if (StaticUtils.isAPermutedVersion(faId.replaceAll("/", "_"), struct.replaceAll("/", "_"))){
//                  found = true;
//                  break;
//                }
//              }
//              if (!found) System.out.println("The positive list does not contain the structure "+className+analyte+" "+faId+"! But we identified it in "+fileName);
//            }else{
//              String analyte = set.getNameStringWithoutRt();
//              if (className.equalsIgnoreCase("TG") && analyte.startsWith("d")) analyte = analyte.substring(1);
//              if (posOfClass==null || (!posOfClass.containsKey(analyte))) System.out.println("The positive list does not contain "+className+" "+msn.getNameStringWithoutRt()+"! But we identified it in "+fileName);
//            }
//          }
//        }
//      }
//    }

      
    }catch (Exception ex){
      ex.printStackTrace();
    }finally{
      if (out!=null){
        try {
          out.close();
        }
        catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
    }   
  }
  
  private Workbook createDetailsExcelFile(String sheetName, int highestNumberOfFiles, int decimalPlacesForMz, double highestMzStdevAllowed,
      double highestRtStdevAllowed, double highestAbundanceCVAllowed, LinkedHashMap<String,Vector<String>> ldaFilesIonMode,
      Hashtable<String,String> appliedMzCorrection, Hashtable<String,LinkedHashMap<String,String>> mods, Hashtable<String,Range> preElutionRanges,
      LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> positiveList, Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> theoreticalMasses,
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>> results, boolean isExperiment1, boolean isExperiment3) throws Exception{
    //this writes the output in the Excel format
    Workbook resultWorkbook = new XSSFWorkbook();
    CellStyle headerStyle = getHeaderStyle(resultWorkbook);
    CellStyle checkStyle = getNotFoundStyle(resultWorkbook);
    CellStyle centerStyle = getCenterStyle(resultWorkbook);
    CellStyle rightIndentStyle = getRightIndentStyle(resultWorkbook);
    Sheet sheet = resultWorkbook.createSheet(sheetName);

    
    int speciesPos = 0;
    String speciesHeader = "Species";
    int longestSpecies = 0;
    int structurePos = 1;
    String structureHeader = "Structure";
    int longestStructure = 0;
    int ionModePos = 2;
    String ionModeHeader = "Ion Mode";
    int longestIonMode = 0;      
    int identificationPos = 3;
    String identificationHeader = "Assigned Structures";
    int longestIdentification = 0;
    int detectionFrequencyPos = 3+highestNumberOfFiles;
    String detectionFrequencyHeader = "Frequency";
    int longestDetectionFrequency = 0;
    int detectionMolFrequencyPos = 4+highestNumberOfFiles;
    String detectionMolFrequencyHeader = "Mol.-Frequency";
    int longestDetectionMolFrequency = 0;
    int formulaPos = 5+highestNumberOfFiles;
    String formulaHeader = "Formula";
    int longestFormula = 0;
    int modPos = 6+highestNumberOfFiles;
    String modHeader = "Adduct";
    int longestMod = 0;
    int theoreticalMzPos = 7+highestNumberOfFiles;
    String theoreticalMzHeader = "Target m/z";
    int longestTheoreticalMz = 0;
    int averageMeasuredMzPos = 8+highestNumberOfFiles;
    String measuredMzHeader = "Measured m/z";
    String averageMeasuredMzHeader = "Avg";
    int longestAverageMeasuredMz = 0;
    int stdevMeasuredMzPos = 9+highestNumberOfFiles;
    String stdevMeasuredMzHeader = "Stdev";
    int longestStdevMeasuredMz = 0;
    int measuredMzValuesPos = 10+highestNumberOfFiles;
    String measuredMzValuesHeader = "Measured m/z values";
    int longestMeasuredMzValues = 0;
    int retentionTimeAveragePos = 10+highestNumberOfFiles*2;
    String retentionTimeHeader = "Retention time";
    String retentionTimeAverageHeader = "Avg";
    int longestRetentionTimeAverage = 0;
    int retentionTimesDeviationPos = 11+highestNumberOfFiles*2;
    String retentionTimesDeviationHeader = "Stdev";
    int longestRetentionTimesDeviation = 0;
    int retentionTimesMeasuredPos = 12+highestNumberOfFiles*2;
    String retentionTimesMeasuredHeader = "Retention times measured";
    int longestRetentionTimesMeasured = 0;
    int commentPos = 12+highestNumberOfFiles*3;
    String commentHeader = "Comment";
    int longestComment = 0;
    int abundancePos = 13+highestNumberOfFiles*3;
    String abundanceHeader = "Abundance";
    String abundanceAverageHeader = "Avg";
    int longestAbundance = 0;
    int abundanceStdevPos = 14+highestNumberOfFiles*3;
    String abundanceStdevHeader = "Stdev";
    int longestAbundanceStdev = 0;
    int abundanceCVPos = 15+highestNumberOfFiles*3;
    String abundanceCVHeader = "CV%";
    int longestAbundanceCV = 0;
    
    int rowCount = 0;
    Row row = null;
    Cell cell = null;
    for (String polarity : ldaFilesIonMode.keySet()){
      int fileNumber = 0;
      for (String filePath : ldaFilesIonMode.get(polarity)){
        String fileName = filePath.substring(filePath.lastIndexOf("\\")+1);
        String massShiftString = "";
        if (!appliedMzCorrection.get(fileName).equalsIgnoreCase("0")) massShiftString = " (massShift="+appliedMzCorrection.get(fileName)+")";
        fileNumber++;
        row = sheet.createRow(rowCount);          
        this.createCell(row, centerStyle, ionModePos, polarity);
        this.createCell(row, null, ionModePos+fileNumber, fileNumber+" = "+fileName+massShiftString);         
        this.createCell(row, centerStyle, measuredMzValuesPos-1, polarity);
        this.createCell(row, null, measuredMzValuesPos+(fileNumber-1), fileNumber+" = "+fileName+massShiftString);
        this.createCell(row, centerStyle, retentionTimesMeasuredPos-1, polarity);
        this.createCell(row, null, retentionTimesMeasuredPos+(fileNumber-1), fileNumber+" = "+fileName+massShiftString);
        
        rowCount++;
      }
    }
    rowCount++;

    int startColumnForRawAreas = abundanceCVPos+1;
    
    row = sheet.createRow(rowCount);
    this.createCell(row, headerStyle, identificationPos, identificationHeader);
    sheet.addMergedRegion(new CellRangeAddress(rowCount,rowCount,identificationPos,identificationPos+highestNumberOfFiles-1));
    this.createCell(row, headerStyle, measuredMzValuesPos, measuredMzValuesHeader);
    sheet.addMergedRegion(new CellRangeAddress(rowCount,rowCount,measuredMzValuesPos,measuredMzValuesPos+highestNumberOfFiles-1));
    this.createCell(row, headerStyle, retentionTimesMeasuredPos, retentionTimesMeasuredHeader);
    sheet.addMergedRegion(new CellRangeAddress(rowCount,rowCount,retentionTimesMeasuredPos,retentionTimesMeasuredPos+highestNumberOfFiles-1));
    this.createCell(row, headerStyle, averageMeasuredMzPos, measuredMzHeader);
    sheet.addMergedRegion(new CellRangeAddress(rowCount,rowCount,averageMeasuredMzPos,averageMeasuredMzPos+1));
    this.createCell(row, headerStyle, retentionTimeAveragePos, retentionTimeHeader);
    sheet.addMergedRegion(new CellRangeAddress(rowCount,rowCount,retentionTimeAveragePos,retentionTimeAveragePos+1));
    this.createCell(row, headerStyle, abundancePos, abundanceHeader);
    sheet.addMergedRegion(new CellRangeAddress(rowCount,rowCount,abundancePos,abundancePos+2));
    rowCount++;
    row = sheet.createRow(rowCount);
    
    this.createCell(row, headerStyle, speciesPos, speciesHeader);
    this.createCell(row, headerStyle, structurePos, structureHeader);
    this.createCell(row, headerStyle, ionModePos, ionModeHeader);
    for (int i=0;i!=highestNumberOfFiles;i++){
      this.createCell(row, headerStyle, identificationPos+i, String.valueOf(i+1));
      this.createCell(row, headerStyle, measuredMzValuesPos+i, String.valueOf(i+1));
      this.createCell(row, headerStyle, retentionTimesMeasuredPos+i, String.valueOf(i+1));
      this.createCell(row, headerStyle, startColumnForRawAreas+i, String.valueOf(i+1));
    }
    this.createCell(row, headerStyle, detectionFrequencyPos, detectionFrequencyHeader);
    this.createCell(row, headerStyle, detectionMolFrequencyPos, detectionMolFrequencyHeader);
    this.createCell(row, headerStyle, formulaPos, formulaHeader);
    this.createCell(row, headerStyle, modPos, modHeader);
    this.createCell(row, headerStyle, theoreticalMzPos, theoreticalMzHeader);
    this.createCell(row, headerStyle, averageMeasuredMzPos, averageMeasuredMzHeader);
    this.createCell(row, headerStyle, stdevMeasuredMzPos, stdevMeasuredMzHeader);
    this.createCell(row, headerStyle, retentionTimeAveragePos, retentionTimeAverageHeader);
    this.createCell(row, headerStyle, retentionTimesDeviationPos, retentionTimesDeviationHeader);
    this.createCell(row, headerStyle, commentPos, commentHeader);
    this.createCell(row, headerStyle, abundancePos, abundanceAverageHeader);
    this.createCell(row, headerStyle, abundanceStdevPos, abundanceStdevHeader);
    this.createCell(row, headerStyle, abundanceCVPos, abundanceCVHeader);
    
    int columnCount = startColumnForRawAreas;
    
    for (String className:positiveList.keySet()) {
      if (!mods.containsKey(className)) continue;
      LinkedHashMap<String,String> modsOfClass= mods.get(className);
      LinkedHashMap<String,LinkedHashMap<String,String>> posAnalytes = positiveList.get(className);
      for (String molName : posAnalytes.keySet()){
        LinkedHashMap<String,String> posStructures = posAnalytes.get(molName);
        Vector<LDAIdentificationDetailsVO> compareVOs = new Vector<LDAIdentificationDetailsVO>();
        if (posStructures.size()>0){
          //if there are no structures present, this hash prevents to add repeatedly the same information
          Hashtable<String,String> noStructuresMod = new  Hashtable<String,String>();
          for (String posStructure:posStructures.keySet()){
            for (LDAIdentificationDetailsVO identVO : extractInfoOfAvailableStructures(sheetName,className,molName,posStructure,posStructures,modsOfClass,ldaFilesIonMode,results,theoreticalMasses,decimalPlacesForMz,preElutionRanges,
                isExperiment1,isExperiment3)){
              if (identVO.containsAnyFileStructures){
                if (identVO.structureFound) compareVOs.add(identVO);
              }else{
                if (!noStructuresMod.containsKey(identVO.originalMod)){
                  compareVOs.add(identVO);
                  noStructuresMod.put(identVO.originalMod, identVO.originalMod);
                }
              }
            }
          }  
        }else{
          compareVOs.addAll(extractInfoOfAvailableStructures(sheetName,className,molName,null,null,modsOfClass,ldaFilesIonMode,results,theoreticalMasses,decimalPlacesForMz,preElutionRanges,isExperiment1,
              isExperiment3));
        }
        for (LDAIdentificationDetailsVO identVO : compareVOs){
          rowCount++;
          row = sheet.createRow(rowCount);
          CellStyle cellStyle = null;
          String comment = null;
          if ((Double.isFinite(identVO.mzStdev) && identVO.mzStdev>highestMzStdevAllowed) ||
              //(Double.isFinite(identVO.rtStdev) && identVO.rtStdev>highestRtStdevAllowed)||
              (Double.isFinite(identVO.areaCV)) && identVO.areaCV>highestAbundanceCVAllowed){
            cellStyle = checkStyle;
            if ((Double.isFinite(identVO.mzStdev) && identVO.mzStdev>highestMzStdevAllowed)){
              if (comment==null) comment = "";
              else comment += ", ";
              comment += "m/z";
            }
            if ( (Double.isFinite(identVO.areaCV)) && identVO.areaCV>highestAbundanceCVAllowed){
              if (comment==null) comment = "";
              else comment += ", ";
              comment += "area";
            }
          }
          for (Double stdev : identVO.rtStdev){
            if (Double.isFinite(stdev) && stdev>highestRtStdevAllowed){
              cellStyle = checkStyle;
              if (comment==null) comment = "";
              else comment += ", ";
              comment += "RT";
            }
          }
          
          String species = className+" "+molName;
          if (species.length()>longestSpecies) longestSpecies = species.length();
          this.createCell(row, cellStyle, speciesPos,species);
          
          String structure = "";
          if (identVO.structure!=null) structure = className+" "+identVO.structure;
          else structure = className+" "+molName;
          if (structure.length()>longestStructure) longestStructure = structure.length();
          this.createCell(row, cellStyle, structurePos, structure);
          
          if (identVO.ionMode.length()>longestMod) longestIonMode = identVO.ionMode.length();
          this.createCell(row, centerStyle, ionModePos, identVO.ionMode);
          
          //if (className.equalsIgnoreCase("PC") && identVO.originalMod.equalsIgnoreCase("Na")) System.out.println(className+" "+molName+" "+structure);            
          String[] idents = identVO.identification.split("\\|");
          //
          if (idents.length!=5 && !(isExperiment1 && sheetName.equalsIgnoreCase("QTRAP 6500") && idents.length==3)){
            System.out.println("1111111111111111111");
            System.out.println(className+molName+" "+structure+": "+idents.length+"    "+identVO.ionMode);
            System.out.println("!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!");
          }
          for (int i=0;i!=idents.length;i++){
            if (idents[i].length()>longestIdentification) longestIdentification = idents[i].length();
            if (!idents[i].equalsIgnoreCase("-"))this.createCell(row, cellStyle, identificationPos+i, idents[i]);              
          }
          if (String.valueOf(identVO.occurence).length()+2>longestDetectionFrequency) longestDetectionFrequency = String.valueOf(identVO.occurence).length()+2;
          this.createCell(row, rightIndentStyle, detectionFrequencyPos, String.valueOf(identVO.occurence));
          
          if (identVO.structOccurrence>0){
            if (String.valueOf(identVO.structOccurrence).length()+2>longestDetectionMolFrequency) longestDetectionMolFrequency = String.valueOf(identVO.structOccurrence).length()+2;
            this.createCell(row, rightIndentStyle, detectionMolFrequencyPos, String.valueOf(identVO.structOccurrence));              
          }
          
          if (identVO.formula.length()>longestFormula) longestFormula = identVO.formula.length();
          this.createCell(row, cellStyle, formulaPos, identVO.formula);
          
          if (identVO.adduct.length()>longestMod) longestMod = identVO.adduct.length();
          this.createCell(row, cellStyle, modPos, identVO.adduct);
          
          String theorMzString = Calculator.FormatNumberToString(identVO.theorMz, decimalPlacesForMz);
          if (theorMzString.length()+2>longestTheoreticalMz) longestTheoreticalMz = theorMzString.length()+2;
          this.createCell(row, rightIndentStyle, theoreticalMzPos, theorMzString);
          
          String measuredMzString = Calculator.FormatNumberToString(identVO.measuredMz, decimalPlacesForMz);
          if (measuredMzString.length()+2>longestAverageMeasuredMz) longestAverageMeasuredMz = measuredMzString.length()+2;
          this.createCell(row, rightIndentStyle, averageMeasuredMzPos, measuredMzString);
          
          String stdevMeasuredMzString = "";
          if (Double.isFinite(identVO.mzStdev)) stdevMeasuredMzString = Calculator.FormatNumberToString(identVO.mzStdev, decimalPlacesForMz+1);
          else stdevMeasuredMzString = "NaN";              
          if (stdevMeasuredMzString.length()+2>longestStdevMeasuredMz) longestStdevMeasuredMz = stdevMeasuredMzString.length()+2;
          this.createCell(row, rightIndentStyle, stdevMeasuredMzPos, stdevMeasuredMzString);
          
          ////if (appliedMzCorrection.get(identVO.ionMode).length()>longestAppliedMzCorrection) longestAppliedMzCorrection = appliedMzCorrection.get(identVO.ionMode).length();
          ////this.createCell(row, cellStyle, appliedMzCorrectionPos, appliedMzCorrection.get(identVO.ionMode));
          
          String[] mzs = identVO.mzs.split("\\|");
          if (idents.length!=5 && !(isExperiment1 && sheetName.equalsIgnoreCase("QTRAP 6500") && idents.length==3)){
            System.out.println("222222222222222222");
            System.out.println("!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!");
          }
          for (int i=0;i!=mzs.length;i++){
            if (mzs[i].length()+2>longestMeasuredMzValues) longestMeasuredMzValues = mzs[i].length()+2;
            if (!mzs[i].equalsIgnoreCase("-"))this.createCell(row, rightIndentStyle, measuredMzValuesPos+i, mzs[i]);              
          }

//          if (identVO.mzs.length()>longestMeasuredMzValues) longestMeasuredMzValues = identVO.mzs.length();
//          this.createCell(row, cellStyle, measuredMzValuesPos, identVO.mzs);
          String rtString = "";
          String stdevRtString = "";
          for (int i=0; i!=identVO.rt.size(); i++){
            if (i>0){
              rtString += "/";
              stdevRtString += "/";
            }
            rtString += identVO.rt.get(i);
            if (Double.isFinite(identVO.rtStdev.get(i))) stdevRtString += Calculator.FormatNumberToString(identVO.rtStdev.get(i), 2);
            else stdevRtString += "NaN";              
          }
          
          if (rtString.length()+2>longestRetentionTimeAverage) longestRetentionTimeAverage= rtString.length()+2;
          this.createCell(row, rightIndentStyle, retentionTimeAveragePos, rtString);

          if (stdevRtString.length()+2>longestRetentionTimesDeviation) longestRetentionTimesDeviation = stdevRtString.length()+2;
          this.createCell(row, rightIndentStyle, retentionTimesDeviationPos, stdevRtString);

//          if (identVO.rts.length()>longestRetentionTimesMeasured) longestRetentionTimesMeasured= identVO.rts.length();
//          this.createCell(row, cellStyle, retentionTimesMeasuredPos, identVO.rts);
          
          String[] rts = identVO.rts.split("\\|");
          if (idents.length!=5 && !(isExperiment1 && sheetName.equalsIgnoreCase("QTRAP 6500") && idents.length==3)){
            System.out.println("333333333333333333");
            System.out.println("!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!");
          }
          for (int i=0;i!=rts.length;i++){
            if (rts[i].length()+2>longestRetentionTimesMeasured) longestRetentionTimesMeasured = rts[i].length()+2;
            if (!rts[i].equalsIgnoreCase("-"))this.createCell(row, rightIndentStyle, retentionTimesMeasuredPos+i, rts[i]);              
          }

          if (comment!=null){
            if (comment.length()+2>longestComment) longestComment = comment.length()+2;
            this.createCell(row, cellStyle, commentPos, comment);              
          }
          String abundanceString = String.valueOf(Math.round(identVO.areaMean));
          if (abundanceString.length()+2>longestAbundance) longestAbundance = abundanceString.length()+2;
          this.createCell(row, rightIndentStyle, abundancePos, abundanceString);
          
          String stdevAreaString = "";
          String cvAreaString = "";
          if (Double.isFinite(identVO.areaStdev)){
            stdevAreaString = String.valueOf(Math.round(identVO.areaStdev));
            cvAreaString = Calculator.FormatNumberToString(identVO.areaCV, 1);
          }else{
            stdevAreaString = "NaN";
            cvAreaString = "NaN";              
          }
          if (stdevAreaString.length()+2>longestAbundanceStdev) longestAbundanceStdev = stdevAreaString.length()+2;
          this.createCell(row, rightIndentStyle, abundanceStdevPos, stdevAreaString);
          if (cvAreaString.length()+2>longestAbundanceCV) longestAbundanceCV = cvAreaString.length()+2;
          this.createCell(row, rightIndentStyle, abundanceCVPos, cvAreaString);
          
          columnCount = startColumnForRawAreas;
          for (int i=0;i!=highestNumberOfFiles;i++){
            String fullPath = ldaFilesIonMode.get(identVO.ionMode).get(i);
            String fileName = fullPath.substring(fullPath.lastIndexOf("\\")+1);
            if (identVO.measuredAreas.containsKey(fileName))
            this.createCell(row, cellStyle, columnCount, String.valueOf(identVO.measuredAreas.get(fileName)));
            columnCount++;
          }
        }
      }
    }
         
    setColumnWidth(sheet, speciesPos, speciesHeader, longestSpecies);
    setColumnWidth(sheet, structurePos, structureHeader, longestStructure);
    setColumnWidth(sheet, ionModePos, ionModeHeader, longestIonMode);
    //setColumnWidth(sheet, identificationPos, identificationHeader, longestIdentification);
    int columnWidth =  (longestIdentification+1)*256;
    for (int i=0; i!=highestNumberOfFiles; i++){
      sheet.setColumnWidth(identificationPos+i,columnWidth);
    }
    setColumnWidth(sheet, detectionFrequencyPos, detectionFrequencyHeader, longestDetectionFrequency);
    setColumnWidth(sheet, detectionMolFrequencyPos, detectionMolFrequencyHeader, longestDetectionMolFrequency);
    setColumnWidth(sheet, formulaPos, formulaHeader, longestFormula);
    setColumnWidth(sheet, modPos, modHeader, longestMod);
    setColumnWidth(sheet, theoreticalMzPos, theoreticalMzHeader, longestTheoreticalMz);
    setColumnWidth(sheet, averageMeasuredMzPos, averageMeasuredMzHeader, longestAverageMeasuredMz);
    setColumnWidth(sheet, stdevMeasuredMzPos, stdevMeasuredMzHeader, longestStdevMeasuredMz);
    //setColumnWidth(sheet, measuredMzValuesPos, measuredMzValuesHeader, longestMeasuredMzValues);
    columnWidth =  (longestMeasuredMzValues+1)*256;
    for (int i=0; i!=highestNumberOfFiles; i++){
      sheet.setColumnWidth(measuredMzValuesPos+i,columnWidth);
    }
    setColumnWidth(sheet, retentionTimeAveragePos, retentionTimesDeviationHeader, longestRetentionTimeAverage);
    setColumnWidth(sheet, retentionTimesDeviationPos, retentionTimesDeviationHeader, longestRetentionTimesDeviation);
    //setColumnWidth(sheet, retentionTimesMeasuredPos, retentionTimesMeasuredHeader, longestRetentionTimesMeasured);   
    columnWidth =  (longestRetentionTimesMeasured+1)*256;
    for (int i=0; i!=highestNumberOfFiles; i++){
      sheet.setColumnWidth(retentionTimesMeasuredPos+i,columnWidth);
    }
    setColumnWidth(sheet, commentPos, commentHeader, longestComment);
    setColumnWidth(sheet, abundancePos, abundanceAverageHeader, longestAbundance);
    setColumnWidth(sheet, abundanceStdevPos, abundanceStdevHeader, longestAbundanceStdev);
    setColumnWidth(sheet, abundanceCVPos, abundanceCVHeader, longestAbundanceCV);
    
    return resultWorkbook;
  }

  private LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> readPositiveList(String file) throws Exception{
    InputStream myxls = new FileInputStream(file);
    Workbook workbook = null;
    if (file.endsWith(".xlsx")) workbook = new XSSFWorkbook(myxls);
    else if (file.endsWith(".xls")) workbook = new HSSFWorkbook(myxls);
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> results = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>>();
    for (int sheetNumber = 0; sheetNumber!=workbook.getNumberOfSheets(); sheetNumber++){
      Sheet sheet = workbook.getSheetAt(sheetNumber);
      String className = sheet.getSheetName();
      LinkedHashMap<String,LinkedHashMap<String,String>> classResults = new LinkedHashMap<String,LinkedHashMap<String,String>>();
      Row row = sheet.getRow(0);
      Vector<Integer> speciesColumn = new Vector<Integer>();
      int structureColumn = -1;
      for (int j=0; j!=row.getLastCellNum(); j++){
        Cell cell = row.getCell(j);
        if (cell==null || cell.getCellType()==Cell.CELL_TYPE_BLANK) continue;
        String contents =  cell.getStringCellValue().trim();
        if (contents.startsWith("MS") && !contents.equalsIgnoreCase("MS/MS"))
          speciesColumn.add(j);
        else if (contents.equalsIgnoreCase("MS/MS")) structureColumn = j;
      }
      if (speciesColumn.size()==0)continue;
      String currentMS1Species = null;
      for (int i=1; i!=(sheet.getLastRowNum()+1); i++){
        row = sheet.getRow(i);
        if (row==null) continue;
        for (Integer j:speciesColumn){
          Cell cell = row.getCell(j);
          if (cell==null || cell.getCellType()==Cell.CELL_TYPE_BLANK) continue;
          String species =  cell.getStringCellValue().trim();
          if (species==null || species.length()==0) continue;
          currentMS1Species = species;
        }
        if (currentMS1Species==null) continue;
        LinkedHashMap<String,String> structures = new LinkedHashMap<String,String>();
        if (classResults.containsKey(currentMS1Species)) structures = classResults.get(currentMS1Species);
        if (structureColumn!=-1){
          Cell cell = row.getCell(structureColumn);
          if (cell!=null && cell.getCellType()!=Cell.CELL_TYPE_BLANK){
            String structure = cell.getStringCellValue().trim();
            if (structure!=null && structure.length()>0){
              if (structure.indexOf(",")!=-1 || structure.indexOf(";")!=-1){
                StringTokenizer tokenizer = new StringTokenizer(structure,",; ");
                while (tokenizer.hasMoreTokens()){
                  String oneStruct = tokenizer.nextToken().trim();
                  structures.put(oneStruct,oneStruct);
                }
              }else{
                structures.put(structure, structure);
              }
            }
          }
        }
        classResults.put(currentMS1Species, structures);
      }
      
      if (classResults.size()>0) results.put(className, classResults);
    }
    
    myxls.close();
    return results;
  }
  
  private void generateDetectedSpeciesListForTG4000QTRAP(){
    String outputFile = "D:\\BiologicalExperiment\\QTRAP\\SpeciesDetectable_MSMS_QTRAP_new.xlsx";
    Vector<String> evaluationFiles = new Vector<String>();
    evaluationFiles.add("D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\positive\\Data20151002_QTrap_Liver-002_QTrap_Liver1-1_pos_LB10_comp.xlsx");
    evaluationFiles.add("D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\positive\\Data20151002_QTrap_Liver-003_QTrap_Liver1-1_pos_LB10_comp.xlsx");
    evaluationFiles.add("D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\positive\\Data20151002_QTrap_Liver-004_QTrap_Liver1-1_pos_LB10_comp.xlsx");
    evaluationFiles.add("D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\positive\\Data20151002_QTrap_Liver-005_QTrap_Liver1-1_pos_LB10_comp.xlsx");
    evaluationFiles.add("D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\positive\\Data20151002_QTrap_Liver-006_QTrap_Liver1-1_pos_LB10_comp.xlsx");
    BufferedOutputStream out = null;
    try {

      Vector<Hashtable<String,Hashtable<String,String>>> foundSpecies = this.parseEvaluationResults(evaluationFiles);
      Hashtable<String,String> coveredMS1Species = new  Hashtable<String,String>();
      LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> possibleSpecies = FoundBiologicalSpecies.getTGSpecies4000QTRAP();
      
      Workbook resultWorkbook = new XSSFWorkbook();
      Sheet sheet = resultWorkbook.createSheet("TG");
      int rowCount = 0;
      Row row;
      for (String species : possibleSpecies.keySet()){
        LinkedHashMap<String,ReferenceInfoVO> structures = possibleSpecies.get(species);
        Hashtable<String,String> coveredMS2Species = new  Hashtable<String,String>();
        boolean speciesWritten = false;
        for (String structure : structures.keySet()){
          //now figure out if the structure was detected
          for (Hashtable<String,Hashtable<String,String>> oneFile : foundSpecies){
            if (!oneFile.containsKey(species))continue;
            coveredMS1Species.put(species, species);
            Hashtable<String,String> foundStructures = oneFile.get(species);
            boolean found = false;
            if (foundStructures.containsKey(structure)) found = true;
            else {
              for (String foundStructure : foundStructures.keySet()){
                if (StaticUtils.isAPermutedVersion(structure, foundStructure,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
                  found = true;
                  break;
                }
              }
            }
            if (found){
              row = sheet.createRow(rowCount);
              if (!speciesWritten){
                Cell cell = row.createCell(0);
                cell.setCellValue(species);
                speciesWritten=true;
              }
              Cell cell = row.createCell(1);
              cell.setCellValue(structure);
              rowCount++;
              coveredMS2Species.put(structure, structure);
              break;
            }
          }
        }
        //cross check if all the detected molecular species are listed in the FoundBiologicalSpecies
        for (int i=0; i!=foundSpecies.size(); i++){
          Hashtable<String,Hashtable<String,String>> oneFile = foundSpecies.get(i);
          if (!oneFile.containsKey(species))continue;
          Hashtable<String,String> foundStructures = oneFile.get(species);
          for (String structure : foundStructures.keySet()){
            boolean found = false;
            if (coveredMS2Species.containsKey(structure)) found = true;
            else {
              for (String foundStructure : coveredMS2Species.keySet()){
                if (StaticUtils.isAPermutedVersion(structure, foundStructure,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
                  found = true;
                  break;
                }
              }
            }
            if (!found){
              System.out.println("The species \""+species+" "+structure+"\" of file "+i+" is not in the FoundBiologicalSpecies!");
            }
          }
        }
      }
      
      //cross check if all the detected species are listed in the FoundBiologicalSpecies
      for (int i=0; i!=foundSpecies.size(); i++){
        Hashtable<String,Hashtable<String,String>> oneFile = foundSpecies.get(i);
        for (String species : oneFile.keySet()){
          Hashtable<String,String> molSpecies = oneFile.get(species);
          if (!coveredMS1Species.containsKey(species)){
            System.out.println("The species \""+species+"\" of file "+i+" is not in the FoundBiologicalSpecies!");
          }
        }
      }
      out = new BufferedOutputStream(new FileOutputStream(outputFile));
      resultWorkbook.close();
      resultWorkbook.write(out);
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } finally {
      try {
        out.close();
      }
      catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
  }
  
  private Vector<Hashtable<String,Hashtable<String,String>>> parseEvaluationResults(Vector<String> filenames) throws Exception{
    Vector<Hashtable<String,Hashtable<String,String>>> allResults = new Vector<Hashtable<String,Hashtable<String,String>>>();
    for (String filename : filenames){
      allResults.add(this.parseEvaluationResult(filename));
    }
    return allResults;
  }
  
  private Hashtable<String,Hashtable<String,String>> parseEvaluationResult(String filename) throws Exception{
    Hashtable<String,Hashtable<String,String>> results = new Hashtable<String,Hashtable<String,String>>();
    InputStream myxls = new FileInputStream(filename);
    Workbook workbook = null;
    workbook = new XSSFWorkbook(myxls);
    for (int i=0; i!=workbook.getNumberOfSheets(); i++){
      String className = workbook.getSheetName(i);
      if (!className.equalsIgnoreCase("TG")) continue;
      Sheet sheet = workbook.getSheetAt(i);
      
      Row row = sheet.getRow(0);
      int speciesColumn = -1;
      int structureColumn = -1;
      int ldaIdentificationCode = -1;
      for (int j=0; j!=row.getLastCellNum(); j++){
        Cell cell = row.getCell(j);
        String contents =  cell.getStringCellValue();
        if (contents.equalsIgnoreCase("LDA")){
          if (speciesColumn<0) speciesColumn = j;
          else structureColumn = j;
        } else if (contents.equalsIgnoreCase("LDA-Code")){
          ldaIdentificationCode = j;
        }
      }
      if (speciesColumn==-1 || structureColumn==-1 || ldaIdentificationCode==-1){
        System.out.println("There is something wrong with the file: "+filename);
        continue;
      }
      
      String currentSpecies = "";
      for (int j=1; j!=(sheet.getLastRowNum()+1);j++){
        row = sheet.getRow(j);
        if (row==null) continue;
        Cell cell = row.getCell(speciesColumn);
        if (cell!=null && cell.getCellType()!=Cell.CELL_TYPE_BLANK){
          String speciesName = cell.getStringCellValue();
          if (!speciesName.equalsIgnoreCase("not reported")) currentSpecies = new String(speciesName);
        }
        
        Cell structureCell = row.getCell(structureColumn);
        Cell identCodeCell = row.getCell(ldaIdentificationCode);
        if (structureCell==null || structureCell.getCellType()==Cell.CELL_TYPE_BLANK || identCodeCell==null || identCodeCell.getCellType()==Cell.CELL_TYPE_BLANK)continue;
        
        String structure = structureCell.getStringCellValue();
        int code = (new Double(identCodeCell.getNumericCellValue())).intValue();
        if (code>0){
          Hashtable<String,String> resultOfOneSpecies = new Hashtable<String,String>();
          if (results.containsKey(currentSpecies)) resultOfOneSpecies = results.get(currentSpecies);
          resultOfOneSpecies.put(structure, structure);
          results.put(currentSpecies, resultOfOneSpecies);
        }
      }
    }
    myxls.close();
    return results;
  }
  
  private void createCell(Row row, CellStyle style, int pos, String value){
    Cell cell = row.createCell(pos);
    cell.setCellValue(value);
    if (cell!=null) cell.setCellStyle(style);
  }
  
  private void setColumnWidth(Sheet sheet, int column, String headerValue, int longestValue){
    int columnWidth = (int)((headerValue.length()*256)*ExcelUtils.BOLD_MULT);
    if ((longestValue+1)*256>columnWidth) columnWidth =  (longestValue+1)*256;
    sheet.setColumnWidth(column,columnWidth); 
  }
  
  @SuppressWarnings("unchecked")
  private Vector<LDAIdentificationDetailsVO> extractInfoOfAvailableStructures(String machineName, String className, String molName, String structure,LinkedHashMap<String,String> posStructures,
      LinkedHashMap<String,String> modsOfClass, LinkedHashMap<String,Vector<String>> ldaFilesIonMode, Hashtable<String,Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>>> results,
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> theoreticalMasses, int mzDecPlaces, Hashtable<String,Range> preElutionRanges, boolean isExperiment1,
      boolean isExperiment3) throws Exception{
    Vector<LDAIdentificationDetailsVO> detailsVOs = new Vector<LDAIdentificationDetailsVO>();
    for (String mod : modsOfClass.keySet()){
      String identification = "";
      String mzs = "";
      String rtString = "";
      String eliteFirstFileNegativeRt = "";
      boolean anyFound = false;
      int occurrence = 0;
      String formula = null;
      boolean containsAnyFileStructures = false;
      boolean structureFound = false;
      Vector<Double> measuredMzs = new Vector<Double>();
      Hashtable<String,Double> measuredAreas = new Hashtable<String,Double>();
      Vector<Vector<Double>> rts = new Vector<Vector<Double>>();
      Hashtable<String,Boolean> structuresFound = new Hashtable<String,Boolean>();
      rts.add(new Vector<Double>());
      Range preElutionRange = null;
      if (preElutionRanges.containsKey(className)){
        preElutionRange = preElutionRanges.get(className);
        rts.add(new Vector<Double>());
      }
      String ionMode = modsOfClass.get(mod);//FoundBiologicalSpecies.getIonModeOfAdduct(mod);
      if (!ldaFilesIonMode.containsKey(ionMode)) continue;
      for (int i=0; i!=ldaFilesIonMode.get(ionMode).size(); i++){
        String filePath = ldaFilesIonMode.get(ionMode).get(i);
        String fileName = filePath.substring(filePath.lastIndexOf("\\")+1);
        Hashtable<String,Hashtable<String,Vector<LipidomicsMSnSet>>> fileResults = results.get(fileName);
        boolean found = false;
        boolean containsStructures = false;
        if (fileResults.containsKey(className) && fileResults.get(className).containsKey(molName)){
          Vector<LipidomicsMSnSet> molResults = fileResults.get(className).get(molName);
          Vector<LipidomicsMSnSet> foundOfMod = new Vector<LipidomicsMSnSet>();
          for (LipidomicsMSnSet msn : molResults){
            if(msn.getModificationName().equalsIgnoreCase(mod)) foundOfMod.add(msn);
          }
          String struct = null;
          if (structure!=null) struct = structure.replaceAll("/", "_");
          if (foundOfMod.size()>0){
            Vector<LipidomicsMSnSet> myResults = new Vector<LipidomicsMSnSet>();
            boolean containsThisFileCorrectStructure = false;
            if (structure!=null){
              //if it contains correct structures, and no FPs, use only the ones with structure information
              // and not the info of additional info without any structure information (only within one file)
              for (LipidomicsMSnSet msn : foundOfMod){
                if (msn.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED){
                  for (Object msnNames : msn.getMSnIdentificationNames()){    
                    String faId = getFaId(msnNames,machineName,className,mod,isExperiment1);
                    //this is to prevent an exclusion, if only FPs are detected
                    boolean matchesACorrectStructure = false;
                    for (String posStruct : posStructures.keySet()){
                      if (StaticUtils.isAPermutedVersion(faId, posStruct.replaceAll("/", "_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
                        matchesACorrectStructure = true;
                        break;
                      }
                    }
                    if (matchesACorrectStructure){
                      containsStructures = true;
                      break;
                    }
                  }
                }
              }
              if (containsStructures){
                containsAnyFileStructures = true;
                for (LipidomicsMSnSet msn : foundOfMod){
                  if (msn.getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED){
                    boolean correctStruct = false;
                    for (Object msnNames : msn.getMSnIdentificationNames()){
                      String faId = getFaId(msnNames,machineName,className,mod,isExperiment1);
                      if (StaticUtils.isAPermutedVersion(faId, struct.replaceAll("/", "_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
                        correctStruct = true;
                        structureFound = true;
                        break;
                      }
                    }
                    if (correctStruct){
                      if (machineName.equalsIgnoreCase("G6550A QTOF") && className.equalsIgnoreCase("PC") &&
                          molName.equalsIgnoreCase("36:0")){
                        float rt = Float.parseFloat(msn.getRt());
                        if (11.0f<rt && rt<13.0f){
                          containsThisFileCorrectStructure = true;
                          myResults.add(msn);
                        }
                        
                      } else if (isExperiment3){
                        double area = getCorrespondingExperiment3Area(fileName,machineName,className,molName,structure,posStructures,msn);
                        if (area>0){
                          containsThisFileCorrectStructure = true;
                          myResults.add(msn);
                        }
                      }else{
                        containsThisFileCorrectStructure = true;
                        myResults.add(msn);
                      }
                    }
                  } else if (isExperiment1 && msn.getStatus()==LipidomicsMSnSet.HEAD_GROUP_DETECTED){
                    myResults.add(msn);
                  }
                }
              }else{
                if (!isExperiment3) myResults.addAll(foundOfMod);
              }
            }else{
              if (machineName.equalsIgnoreCase("Q Exactive") && className.equalsIgnoreCase("P-PC")){
                for (LipidomicsMSnSet msn : foundOfMod){
                  float rt = Float.parseFloat(msn.getRt());
                  if (molName.equalsIgnoreCase("38:4") && 28.7f<rt && rt<29.7f)
                    myResults.add(msn);
                }
              }else
                myResults.addAll(foundOfMod);
            }
            if (myResults.size()>0){
              found = true;
              anyFound = true;
              occurrence++;
              if (containsThisFileCorrectStructure){
                structuresFound.put(fileName, true);
              }
              //now I have the relevent ones -> extract the required information
              double highestArea = 0d;
              String strongestAssignment = "";
              double mzOneFile = 0d;
              double monoIntensityOneFile = 0d;
              double areaOneFile = 0d;
              double areaOneFilePreElution = 0d;
              String rt = "";
              String rtPreElution = "";
              double rtOneFile = 0d;
              double rtPreElutionOneFile = 0d;
              
              if (myResults.size()>1 && preElutionRange==null){
                System.out.println("There is more than one peak for "+className+molName+" "+(structure!=null ? structure:"")+"  "+mod+" at file "+fileName);
              }

              for (int j=0;j!=myResults.size();j++){
                LipidomicsMSnSet msn = myResults.get(j);
                formula = StaticUtils.getFormulaInHillNotation(StaticUtils.categorizeAdduct(msn.getAnalyteFormula()), false);
                for (CgProbe probe : msn.getIsotopicProbes().get(0)){
                  mzOneFile += (double)(probe.Mz*probe.Area);
                  monoIntensityOneFile += probe.Area;
                }
                if (containsStructures){
                  for (Object msnNames : msn.getMSnIdentificationNames()){
                    String faId = getFaId(msnNames,machineName,className,mod,isExperiment1);
                    boolean isPermutedVersion = StaticUtils.isAPermutedVersion(faId, struct.replaceAll("/", "_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
                    if (!isPermutedVersion && !isExperiment1) continue;
                    String nameString = "";
                    if (msnNames instanceof Vector){
                      for (String name : (Vector<String>)msnNames){
                        nameString += name+"|";
                      }
                      nameString = nameString.substring(0,nameString.length()-1);
                    }else{
                      nameString = (String)msnNames;
                    }
                    double area = 0;
                    if (isPermutedVersion){
                      if (isExperiment3){
                        area = getCorrespondingExperiment3Area(fileName,machineName,className,molName,structure,posStructures,msn);
                      }else if (isExperiment1)
                        area = msn.getArea();                        
                      else 
                        area = msn.getArea()*msn.getRelativeIntensity(nameString);
                    }else if (isExperiment1 && msn.getStatus()<LipidomicsMSnSet.FRAGMENTS_DETECTED){
                      if (!isExperiment3) area = msn.getArea();
                    }
                    
                    if (myResults.size()>1){
                      if (preElutionRange!=null && preElutionRange.insideRange(Float.parseFloat(msn.getRt()))){
                        rtPreElutionOneFile += Double.valueOf(msn.getRt())*area;
                        areaOneFilePreElution += area;                        
                      }else{
                        rtOneFile += Double.valueOf(msn.getRt())*area;
                        areaOneFile += area;
                      }
                      
                    }else{
                      areaOneFile += area;
                      if (preElutionRange!=null && preElutionRange.insideRange(Float.parseFloat(msn.getRt())))
                        rtPreElution = msn.getRt();
                      else
                        rt = msn.getRt();
                    }
                    if (area>highestArea){
                      highestArea = area;
                      strongestAssignment = nameString.replaceAll("\\|", ",");
                    }
                  }
                }else{
                  if (isExperiment3) continue;
                  if (myResults.size()>1){
                    if (preElutionRange!=null && preElutionRange.insideRange(Float.parseFloat(msn.getRt()))){
                      rtPreElutionOneFile += Double.valueOf(msn.getRt())*msn.getArea();
                      areaOneFilePreElution += msn.getArea();
                    }else{
                      rtOneFile += Double.valueOf(msn.getRt())*msn.getArea();
                      areaOneFile += msn.getArea();
                    }
                  }else{
                    areaOneFile +=msn.getArea();
                    if (preElutionRange!=null && preElutionRange.insideRange(Float.parseFloat(msn.getRt())))
                      rtPreElution = msn.getRt();
                    else
                      rt = msn.getRt();
                  }                  
                }
//                String correctIdent = "";
//                for (Object msnNames : msn.getMSnIdentificationNames()){
//                  
//                }
              }
              
              mzOneFile = mzOneFile/monoIntensityOneFile;
              if (myResults.size()>1){
                if (rtOneFile!=0d)
                  rt = Calculator.FormatNumberToString(rtOneFile/areaOneFile,2);
                if (rtPreElutionOneFile!=0d)
                  rtPreElution = Calculator.FormatNumberToString(rtPreElutionOneFile/areaOneFilePreElution,2);              
              }
              measuredMzs.add(mzOneFile);
              measuredAreas.put(fileName, areaOneFile+areaOneFilePreElution);
              //System.out.println(className+molName+" "+(structure!=null ? structure:"")+"  "+mod+"  "+rt+")");
              if (preElutionRange!=null){
                if (rtPreElution.length()>0){
                  rts.get(0).add(new Double(rtPreElution));
                  rtString += rtPreElution;
                }
                if (rtPreElution.length()>0&&rt.length()>0)
                  rtString += "/";
                if (rt.length()>0){      
                  //this is a special case only for the first file of the Orbitrap Elite
                  if (!isExperiment1 && machineName.equalsIgnoreCase("Orbitrap Elite") && modsOfClass.get(mod).equalsIgnoreCase("-")&&i==0)
                    eliteFirstFileNegativeRt = rt;
                  else
                    rts.get(1).add(new Double(rt));
                  rtString += rt;
                }
              }else{
                //this is a special case only for the first file of the Orbitrap Elite
                if (!isExperiment1 && machineName.equalsIgnoreCase("Orbitrap Elite") && modsOfClass.get(mod).equalsIgnoreCase("-")&&i==0)
                  eliteFirstFileNegativeRt = rt;
                else
                  rts.get(0).add(new Double(rt));                
                rtString += rt;
              } 
              rtString += "|";
              if (containsStructures){
                identification += className+" "+strongestAssignment+"|";
              }else{
                identification += className+" "+molName+"|";
              }
              mzs += Calculator.FormatNumberToString(mzOneFile, mzDecPlaces)+"|";
              
            }
          }
        }
        if (!found){
          identification += "-|";
          mzs += "-|";
          rtString += "-|";
        }
      }
      if (anyFound){
        LDAIdentificationDetailsVO details = new LDAIdentificationDetailsVO();
        String struct = structure;
        if (structure!=null){
          struct = new String(structure);
          if (className.equalsIgnoreCase("DG") && mod.equalsIgnoreCase("NH4")){
            struct = struct.replaceAll("/", "_");
            if (struct.contains("_-")) struct = struct.substring(0,struct.indexOf("_-"))+struct.substring(struct.indexOf("_-")+2);
            else if (struct.contains("-_")) struct = struct.substring(0,struct.indexOf("-_"))+struct.substring(struct.indexOf("-_")+2);
          }// else if ((machineName.equalsIgnoreCase("Orbitrap Velos Pro CID")||machineName.equalsIgnoreCase("Orbitrap Velos Pro HCD")) &&
            //  className.equalsIgnoreCase("DG") && mod.equalsIgnoreCase("Na") && struct.split("_").length==2){
           // struct = struct += "_-";
          //}
        }
        
        Hashtable<String,Double> usableAreas = new Hashtable<String,Double>();
        if (containsAnyFileStructures){
          details.structure = struct;
          for (int i=0; i!=ldaFilesIonMode.get(ionMode).size(); i++){
            String filePath = ldaFilesIonMode.get(ionMode).get(i);
            String fileName = filePath.substring(filePath.lastIndexOf("\\")+1);
            if (structuresFound.containsKey(fileName) || (isExperiment1 && measuredAreas.containsKey(fileName))) usableAreas.put(fileName, measuredAreas.get(fileName));
          }
        } else{
          details.structure = molName;
          usableAreas = new Hashtable<String,Double>(measuredAreas);
        }
        identification = identification.substring(0,identification.length()-1);
        mzs = mzs.substring(0,mzs.length()-1);
        rtString = rtString.substring(0,rtString.length()-1);
        details.identification = identification;
        details.occurence = occurrence;
        if (className.equalsIgnoreCase("LPC") || className.equalsIgnoreCase("LPE") || className.equalsIgnoreCase("LPS")||
            className.equalsIgnoreCase("SM") || className.equalsIgnoreCase("Cer")){
          details.structOccurrence = occurrence;
        }else{
          details.structOccurrence = structuresFound.size();
        }
        details.formula = formula;
        details.originalMod = mod;
        details.adduct = FoundBiologicalSpecies.getIonisation(mod);
        details.ionMode = modsOfClass.get(mod);
        details.theorMz = theoreticalMasses.get(className).get(molName).get(mod).getAnalyteMass();
        details.measuredMz = Calculator.mean(measuredMzs);
        if (measuredMzs.size()>1) details.mzStdev = Calculator.stddeviation(measuredMzs);
        else details.mzStdev = Double.NaN;
        details.mzs = mzs;
        details.containsAnyFileStructures = containsAnyFileStructures;
        details.structureFound = structureFound;
        detailsVOs.add(details);
        details.rt = new Vector<String>();
        details.rtStdev = new Vector<Double>();
        if (!isExperiment1 && machineName.equalsIgnoreCase("Orbitrap Elite") && modsOfClass.get(mod).equalsIgnoreCase("-") && eliteFirstFileNegativeRt.length()>0){
          if (rts.size()>1 && rts.get(0).size()>0){
            details.rt.add(Calculator.FormatNumberToString(Calculator.mean(rts.get(0)),2));
            details.rtStdev.add(Calculator.stddeviation(rts.get(0)));
            details.rt.add(eliteFirstFileNegativeRt);
            details.rtStdev.add(Double.NaN);
            if (rts.get(1).size()>0){
              details.rt.add(Calculator.FormatNumberToString(Calculator.mean(rts.get(1)),2));
              details.rtStdev.add(Calculator.stddeviation(rts.get(1)));              
            }
          }else{
            details.rt.add(eliteFirstFileNegativeRt);
            details.rtStdev.add(Double.NaN);
            if (rts.get(0).size()>0){
              details.rt.add(Calculator.FormatNumberToString(Calculator.mean(rts.get(0)),2));
              details.rtStdev.add(Calculator.stddeviation(rts.get(0)));
            }            
          }
        }else{
          if (rts.get(0).size()>0){
            details.rt.add(Calculator.FormatNumberToString(Calculator.mean(rts.get(0)),2));
            details.rtStdev.add(Calculator.stddeviation(rts.get(0)));
          }
          if (rts.size()>1&&rts.get(1).size()>0){
            details.rt.add(Calculator.FormatNumberToString(Calculator.mean(rts.get(1)),2));
            details.rtStdev.add(Calculator.stddeviation(rts.get(1)));          
          }
        }
        details.rts = rtString;
        details.areaMean = Calculator.mean(new Vector<Double>(usableAreas.values()));
        details.areaStdev = Calculator.stddeviation(new Vector<Double>(usableAreas.values()));
        if (Double.isFinite(details.areaStdev)) details.areaCV = (details.areaStdev*100d)/details.areaMean;
        details.measuredAreas = usableAreas;
      }
    }
    return detailsVOs;
  }
  
  @SuppressWarnings("unchecked")
  private String getFaId(Object msnNames, String machineName, String lipidClass, String mod, boolean isExperiment1){
    String nameString = "";
    String faId = "";
    if (msnNames instanceof Vector){
      for (String name : (Vector<String>)msnNames){
        nameString += name+"|";
      }
      faId = nameString.substring(0,nameString.indexOf("|"));
      nameString = nameString.substring(0,nameString.length()-1);
    }else{
      nameString = (String)msnNames;
      faId = nameString;
    }
    faId = faId.replaceAll("/", "_");
    if (isExperiment1){
      if (lipidClass.equalsIgnoreCase("DG") && faId.split("_").length==2){
        faId+="_-";
      }      
    }else{
      if ((machineName.equalsIgnoreCase("Orbitrap Velos Pro CID")||machineName.equalsIgnoreCase("Orbitrap Velos Pro HCD")) &&
          ((lipidClass.equalsIgnoreCase("DG") && mod.equalsIgnoreCase("NH4")) ||
           (lipidClass.equalsIgnoreCase("DG") && mod.equalsIgnoreCase("Na") && faId.split("_").length==2))){
        faId+="_-";
      }
    }
    return faId;
  }
  
  private class LDAIdentificationDetailsVO{
    public String structure = null;
    public String identification = null;
    public int occurence = -1;
    public int structOccurrence = -1;
    public String formula = null;
    public String originalMod = null;
    public String adduct = null;
    public String ionMode = null;
    public double theorMz = 0;
    public double measuredMz = 0;
    public double mzStdev = 0;
    public String mzs = null;
    public boolean containsAnyFileStructures = false;
    public boolean structureFound = false;
    public Vector<String> rt = null;
    public Vector<Double> rtStdev = null;
    public String rts = null;
    public double areaMean = 0;
    public double areaStdev = 0;
    public double areaCV = 0;
    public Hashtable<String,Double> measuredAreas;
  }
  
  private void calcMeanAndStdev(){
    Vector<Double> values = new Vector<Double>();
    values.add(19.00d);
    values.add(19.10d);
//    values.add(32.70d);
//    values.add(31.39d);
//    values.add(30.51d);
    System.out.println("Mean:\t"+Calculator.mean(values));
    System.out.println("Stdev:\t"+Calculator.stddeviation(values));
  }
  
  private void crossPlatformComparison(){
    try{
      String detailsFileName = "D:\\BiologicalExperiment\\SupplementaryTable9.xlsx";
      String speciesSequenceFile = "D:\\BiologicalExperiment\\massLists\\positive\\positive.xlsx";
      String resultFile = "D:\\BiologicalExperiment\\CrossPlatformComparison_new.xlsx";
      String novelFile = "D:\\BiologicalExperiment\\NovelSpecies.xlsx";
      
      LinkedHashMap<String,Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>>> resultDetails = readDetailsOfSeveralPlatforms(detailsFileName);
      Hashtable<String,Vector<String>> analyteSequence = (Hashtable<String,Vector<String>>)QuantificationThread.parseQuantExcelFile(speciesSequenceFile,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f,false).get(1);
      Vector<String> classSequence = new Vector<String>();
      classSequence.add("PI");
      classSequence.add("P-PC");
      classSequence.add("P-PE");
      classSequence.add("LPC");
      classSequence.add("LPE");
      classSequence.add("PS");
      classSequence.add("LPS");
      classSequence.add("PC");
      classSequence.add("PE");
      classSequence.add("PG");
      classSequence.add("DG");
      classSequence.add("TG");
      classSequence.add("SM");
      classSequence.add("Cer");     
      writeCrossPlatformComparisonFile(resultFile,novelFile,classSequence,analyteSequence,resultDetails);
    }catch (Exception ex){
      
    }
  }
  
  private LinkedHashMap<String,Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>>> readDetailsOfSeveralPlatforms(String fileName) throws Exception{
    InputStream myxls = null;
    Workbook workbook = null;
    LinkedHashMap<String,Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>>> verifiedResults = new LinkedHashMap<String,Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>>>();
    try{
      myxls = new FileInputStream(fileName);      
      workbook = new XSSFWorkbook(myxls);

      for (int sheetNr=0; sheetNr!=workbook.getNumberOfSheets(); sheetNr++){
        Sheet sheet = workbook.getSheetAt(sheetNr);
        Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>> oneSheet = readDetailsOfOnePlatform(sheet);
        if (oneSheet.size()>0) verifiedResults.put(sheet.getSheetName(), oneSheet);
      }
      
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } finally {
      try {
        if (workbook!=null) workbook.close();
      } catch (IOException e) {e.printStackTrace();}
      try {
        if (myxls!=null) myxls.close();
      } catch (IOException e) {e.printStackTrace();}
    }
    return verifiedResults;
  }
  
  private Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>> readDetailsOfOnePlatform(Sheet sheet) throws Exception{
    Row row = null;
    
    boolean topRowFound = false;
    boolean detailsRowFound = false;
    
    RangeInteger structuresColumns = null;
    RangeInteger measuredMzColumns = null;
    RangeInteger retentionTimeColumns = null;
    
    int speciesColumn = -1;
    int structureColumn = -1;
    int ionModeColumn = -1;
    int frequencyColumn = -1;
    int molFrequencyColumn = -1;
    int formulaColumn = -1;
    int adductColumn = -1;
    int targetMzColumn = -1;
    int measuredMzColumn = -1;
    //int retentionTimeColumn = -1;
    
    boolean writtenOutOnce = false;
    
    Hashtable<String,Range> preElutionRanges = new Hashtable<String,Range>();
    if (sheet.getSheetName().equalsIgnoreCase("G6550A QTOF"))
      preElutionRanges = FoundBiologicalSpecies.getQTOFG6550PreElutionRanges();
    else if (sheet.getSheetName().equalsIgnoreCase("Orbitrap Elite"))
      preElutionRanges = FoundBiologicalSpecies.getOrbitrapElitePreElutionRanges();
 
    
    Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>> results = new Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>>();
    for (int i=0; i!=(sheet.getLastRowNum()+1);i++){
      row = sheet.getRow(i);
      if (row==null) continue;
            
      if (!topRowFound || !detailsRowFound){
        if (!topRowFound){
          for (int j=0; j!=row.getLastCellNum(); j++){
            Cell cell = row.getCell(j);
          
            String contents = "";
            Double numeric = null;
            int cellType = -1;
            if (cell!=null) cellType = cell.getCellType();
            if (cellType==Cell.CELL_TYPE_STRING){
              contents = cell.getStringCellValue().trim();
              try{ 
                if (contents!=null)numeric = new Double(contents.replaceAll(",", "."));
              }catch(NumberFormatException nfx){};
            }else if (cellType==Cell.CELL_TYPE_NUMERIC || cellType==Cell.CELL_TYPE_FORMULA){
              numeric = cell.getNumericCellValue();
              contents = String.valueOf(numeric);
            }
            if (contents.equalsIgnoreCase("Assigned Structures")){
              int endColumn = getEndColumnOfMergedRegion(sheet, i, j);
              if (endColumn>-1) structuresColumns = new RangeInteger(j,endColumn);
            }
            if (contents.equalsIgnoreCase("Measured m/z")){
              int endColumn = getEndColumnOfMergedRegion(sheet, i, j);
              if (endColumn>-1) measuredMzColumns = new RangeInteger(j,endColumn);
            }
            if (contents.startsWith("Retention times measured")){
              int endColumn = getEndColumnOfMergedRegion(sheet, i, j);
              if (endColumn>-1) retentionTimeColumns = new RangeInteger(j,endColumn);
            }
          }
          if (structuresColumns!=null && measuredMzColumns!=null && retentionTimeColumns!=null) topRowFound = true;
        } else {
          for (int j=0; j!=row.getLastCellNum(); j++){
            Cell cell = row.getCell(j);
          
            String contents = "";
            Double numeric = null;
            int cellType = -1;
            if (cell!=null) cellType = cell.getCellType();
            if (cellType==Cell.CELL_TYPE_STRING){
              contents = cell.getStringCellValue().trim();
              try{ 
                if (contents!=null)numeric = new Double(contents.replaceAll(",", "."));
              }catch(NumberFormatException nfx){};
            }else if (cellType==Cell.CELL_TYPE_NUMERIC || cellType==Cell.CELL_TYPE_FORMULA){
              numeric = cell.getNumericCellValue();
              contents = String.valueOf(numeric);
            }
            
            if (contents.equalsIgnoreCase("Species"))
              speciesColumn = j;
            else if (contents.equalsIgnoreCase("Structure"))
              structureColumn = j;
            else if (contents.equalsIgnoreCase("Ion Mode"))
                ionModeColumn = j;
            else if (contents.equalsIgnoreCase("Frequency"))
              frequencyColumn = j;
            else if (contents.equalsIgnoreCase("Mol.-Frequency"))
              molFrequencyColumn = j;
            else if (contents.equalsIgnoreCase("Formula"))
              formulaColumn = j;
            else if (contents.equalsIgnoreCase("Adduct"))
              adductColumn = j;
            else if (contents.equalsIgnoreCase("Target m/z"))
              targetMzColumn = j;
            else if (measuredMzColumns.insideRange(j) && contents.equalsIgnoreCase("Avg"))
              measuredMzColumn = j;
//            else if (retentionTimeColumns.insideRange(j) && contents.equalsIgnoreCase("Avg"))
//              retentionTimeColumn = j;
          }
          if (speciesColumn>-1 && structureColumn>-1 && ionModeColumn>-1 && frequencyColumn>-1
              && molFrequencyColumn>-1 && formulaColumn>-1 && adductColumn>-1 && targetMzColumn>-1
              && measuredMzColumn>-1)//&& retentionTimeColumn>-1)
            detailsRowFound = true;
        }
      }else{
        Cell cell = row.getCell(speciesColumn);
        if (cell==null || cell.getCellType()==Cell.CELL_TYPE_BLANK) continue;
        SingleAdductIdentificationVO identVO = new SingleAdductIdentificationVO();
        identVO.species = getStringOfMandatoryString(cell, sheet.getSheetName(), i);
        cell = row.getCell(structureColumn);
        identVO.molecularSpecies = getStringOfMandatoryString(cell, sheet.getSheetName(), i);
        if (identVO.species.startsWith("DG ")){
          if ((identVO.molecularSpecies.indexOf("_")!=-1 && identVO.molecularSpecies.indexOf("_")==identVO.molecularSpecies.lastIndexOf("_"))||
            (identVO.molecularSpecies.indexOf("/")!=-1 && identVO.molecularSpecies.indexOf("/")==identVO.molecularSpecies.lastIndexOf("/"))){
            identVO.molecularSpecies = identVO.molecularSpecies.replaceAll("/", "_");
            identVO.molecularSpecies += "_-";
          }
        }
        cell = row.getCell(ionModeColumn);
        identVO.ionMode = getStringOfMandatoryString(cell, sheet.getSheetName(), i);
        cell = row.getCell(adductColumn);
        identVO.adduct = getStringOfMandatoryString(cell, sheet.getSheetName(), i);
        for (int k=structuresColumns.getStart(); k!=(structuresColumns.getStop()+1); k++){
          cell = row.getCell(k);
          if (cell==null || cell.getCellType()==Cell.CELL_TYPE_BLANK) continue;
          String assignedStructure = cell.getStringCellValue();
          int msRunNumber = k-structuresColumns.getStart()+1;
          if (assignedStructure!=null && assignedStructure.length()>0){
            assignedStructure = assignedStructure.trim();
            if (assignedStructure.startsWith("DG ") && identVO.adduct.equalsIgnoreCase("+Na+")
                && assignedStructure.indexOf("/")!=-1 && assignedStructure.indexOf("_")!=-1)
              assignedStructure = assignedStructure.replaceAll("/", "_");
            identVO.annotatedStructures.put(msRunNumber, assignedStructure);

          }
        }
        cell = row.getCell(frequencyColumn);
        identVO.frequency = getIntegerCellValue(cell, sheet.getSheetName(), i, true);
        cell = row.getCell(molFrequencyColumn);
        identVO.molFrequency = getIntegerCellValue(cell, sheet.getSheetName(), i, false);
        cell = row.getCell(formulaColumn);
        identVO.formula = getStringOfMandatoryString(cell, sheet.getSheetName(), i);
        cell = row.getCell(targetMzColumn);
        identVO.targetMz = getStringOfMandatoryString(cell, sheet.getSheetName(), i);
        cell = row.getCell(measuredMzColumn);
        identVO.measuredMz = "";
        if (cell.getCellType()==Cell.CELL_TYPE_STRING)
          identVO.measuredMz = getStringOfMandatoryString(cell, sheet.getSheetName(), i);
        else if (cell.getCellType()==Cell.CELL_TYPE_NUMERIC)
          identVO.measuredMz = Calculator.FormatNumberToString(cell.getNumericCellValue(), identVO.targetMz.length()-identVO.targetMz.indexOf(".")-1);
        for (int k=retentionTimeColumns.getStart(); k!=(retentionTimeColumns.getStop()+1); k++){
          cell = row.getCell(k);
          if (cell==null || cell.getCellType()==Cell.CELL_TYPE_BLANK) continue;
          String rtString = null;
          if (cell.getCellType()==Cell.CELL_TYPE_STRING)
            rtString = cell.getStringCellValue();
          else if (cell.getCellType()==Cell.CELL_TYPE_NUMERIC)
            rtString = Calculator.FormatNumberToString(cell.getNumericCellValue(),2);
          else throw new Exception("An RT entry must be numeric or String; the entry at sheet "+sheet.getSheetName()+" row number "+(i+1)+" is not!");
          if (rtString==null && rtString.length()==0) continue;
          int msRunNumber = k-retentionTimeColumns.getStart()+1;
          if (sheet.getSheetName().equalsIgnoreCase("G6550A QTOF") || sheet.getSheetName().equalsIgnoreCase("Orbitrap Elite")){
            Double firstElution = null;
            Double secondElution = null;
            String lipidClass = identVO.species.substring(0,identVO.species.indexOf(" "));
            if (rtString.indexOf("/")!=-1){
              firstElution = new Double(rtString.substring(0,rtString.indexOf("/")));
              secondElution = new Double(rtString.substring(rtString.indexOf("/")+1));
            } else firstElution = new Double(rtString);
            if (secondElution!=null && !preElutionRanges.containsKey(lipidClass))
              throw new Exception("This splitted RT entry is not in the pre-elution ranges "+sheet.getSheetName()+" row number "+(i+1)+"!");
            if (sheet.getSheetName().equalsIgnoreCase("G6550A QTOF")){
              if (preElutionRanges.containsKey(lipidClass)){
                if (preElutionRanges.get(lipidClass).insideRange(firstElution.floatValue())){
                  identVO.firstRetentionTimes.put(msRunNumber, firstElution);
                  if (secondElution!=null) identVO.secondRetentionTimes.put(msRunNumber, secondElution);
                } else {
                  identVO.secondRetentionTimes.put(msRunNumber, firstElution);
                  if (secondElution!=null) throw new Exception("There cannot be a splitted RT if the first RT is not inside the pre-elution ranges: "+sheet.getSheetName()+" row number "+(i+1)+"!");
                }
              }else{
                identVO.secondRetentionTimes.put(msRunNumber, firstElution);
              }
            } else if (sheet.getSheetName().equalsIgnoreCase("Orbitrap Elite")){
              if (identVO.ionMode.equalsIgnoreCase("+") || (identVO.ionMode.equalsIgnoreCase("-") && msRunNumber==1)){
                if (preElutionRanges.containsKey(lipidClass)){
                  if (preElutionRanges.get(lipidClass).insideRange(firstElution.floatValue())){
                    identVO.firstRetentionTimes.put(msRunNumber, firstElution);
                    if (secondElution!=null) identVO.secondRetentionTimes.put(msRunNumber, secondElution);
                  } else {
                    identVO.secondRetentionTimes.put(msRunNumber, firstElution);
                    if (secondElution!=null) throw new Exception("There cannot be a splitted RT if the first RT is not inside the pre-elution ranges: "+sheet.getSheetName()+" row number "+(i+1)+"!");
                  }
                }else{
                  identVO.secondRetentionTimes.put(msRunNumber, firstElution);
                }                
              } else if (identVO.ionMode.equalsIgnoreCase("-") && msRunNumber!=1){
                if (secondElution!=null) throw new Exception("It is not possible that a split occurs for negative ion mode runs 2-4: "+sheet.getSheetName()+" row number "+(i+1)+"!");
                identVO.thirdRetentionTimes.put(msRunNumber, firstElution);
              } else  throw new Exception("This composition is not possible: "+sheet.getSheetName()+" row number "+(i+1)+"!");
            }
          }else{
            identVO.firstRetentionTimes.put(msRunNumber, new Double(rtString));
          }
        }

////        cell = row.getCell(retentionTimeColumn);
////        identVO.retentionTime = "";
      ////        if (cell.getCellType()==Cell.CELL_TYPE_STRING)
      ////          identVO.retentionTime = getStringOfMandatoryString(cell, sheet.getSheetName(), i);
      ////        else if (cell.getCellType()==Cell.CELL_TYPE_NUMERIC)
      ////        identVO.retentionTime = Calculator.FormatNumberToString(cell.getNumericCellValue(),2);
        
        LinkedHashMap<String,Vector<SingleAdductIdentificationVO>> speciesDetections = new LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>();
        if (results.containsKey(identVO.species)) speciesDetections = results.get(identVO.species);
        Vector<SingleAdductIdentificationVO> structureIdents = new Vector<SingleAdductIdentificationVO>();
        if (speciesDetections.containsKey(identVO.molecularSpecies)) structureIdents = speciesDetections.get(identVO.molecularSpecies);
        else if (identVO.species.startsWith("DG ")){
          for (String molSpecies:speciesDetections.keySet()){
            String one = molSpecies.replaceAll("/", "_").substring("DG ".length());
            String two = identVO.molecularSpecies.replaceAll("/", "_").substring("DG ".length());
            if (!StaticUtils.isAPermutedVersion(one, two,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) continue;
            identVO.molecularSpecies = molSpecies;
            structureIdents = speciesDetections.get(identVO.molecularSpecies);
            break;
          }
        }
        structureIdents.add(identVO);
        speciesDetections.put(identVO.molecularSpecies, structureIdents);
        results.put(identVO.species,speciesDetections);
      }
    }
    return results;
  }
  
  private int getEndColumnOfMergedRegion(Sheet sheet, int startRow, int startColumn){
    for (int k=0; k!=sheet.getNumMergedRegions(); k++){
      CellRangeAddress cra =  sheet.getMergedRegion(k);
      if (cra.getFirstRow()!=startRow || cra.getFirstColumn()!=startColumn) continue;
      return cra.getLastColumn();              
    }
    return -1;
  }
  
  private String getStringOfMandatoryString(Cell cell, String sheetName, int rowNr) throws Exception{
    if (cell==null || cell.getCellType()==Cell.CELL_TYPE_BLANK) throw new Exception("There is something wrong in sheet "+sheetName+" at row "+(rowNr+1));
    String value = cell.getStringCellValue();
    if (value==null || value.length()==0) throw new Exception("There is something wrong in sheet "+sheetName+" at row "+(rowNr+1));
    return value;
  }
  
  private Integer getIntegerCellValue(Cell cell, String sheetName, int rowNr, boolean mandatory) throws Exception{
    if (cell==null || cell.getCellType()==Cell.CELL_TYPE_BLANK){
      if (mandatory) throw new Exception("There is something wrong in sheet "+sheetName+" at row "+(rowNr+1));
      else return null;
    }
    String contents = "";
    Integer numeric = null;
    int cellType = -1;
    if (cell!=null) cellType = cell.getCellType();
    if (cellType==Cell.CELL_TYPE_STRING){
      contents = cell.getStringCellValue().trim();
      try{ 
        if (contents!=null) numeric = (int)Math.round(new Double(contents.replaceAll(",", ".")));
      }catch(NumberFormatException nfx){};
    }else if (cellType==Cell.CELL_TYPE_NUMERIC || cellType==Cell.CELL_TYPE_FORMULA){
      numeric = (int)Math.round(cell.getNumericCellValue());
      contents = String.valueOf(numeric);
    }
    return numeric;
  }
  
  private class SingleAdductIdentificationVO{
    
    public SingleAdductIdentificationVO(){
      annotatedStructures = new Hashtable<Integer,String>();
    }
    
    public SingleAdductIdentificationVO(SingleAdductIdentificationVO other){
      this();
      this.species = other.species;
      this.molecularSpecies = other.molecularSpecies;
      this.ionMode = other.ionMode;
      this.frequency = other.frequency;
      this.molFrequency = other.molFrequency;
      this.formula = other.formula;
      this.adduct = other.adduct;
      this.targetMz = other.targetMz;
      this.measuredMz = other.measuredMz;
    }
    
    
    public String species = null;
    public String molecularSpecies = null;
    public String ionMode = null;
    public Hashtable<Integer,String> annotatedStructures = null;
    public Integer frequency = null;
    public Integer molFrequency = null;
    public String formula = null;
    public String adduct = null;
    public String targetMz = null;
    public String measuredMz = null;
    //public String retentionTime = null;
    public Hashtable<Integer,Double> firstRetentionTimes = new Hashtable<Integer,Double>();
    public Hashtable<Integer,Double> secondRetentionTimes = new Hashtable<Integer,Double>();
    public Hashtable<Integer,Double> thirdRetentionTimes = new Hashtable<Integer,Double>();
  }
  
  
  private void writeCrossPlatformComparisonFile(String fileName, String fileNameNovel, Vector<String> classSequence, Hashtable<String,Vector<String>> analyteSequence,
      LinkedHashMap<String,Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>>> resultDetails) throws Exception{
    
    BufferedOutputStream out = null;
    BufferedOutputStream outNovel = null;
    try{
      Workbook resultWorkbook = new XSSFWorkbook();
      Workbook novelSpeciesWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      CellStyle centerStyle = getCenterStyle(resultWorkbook);
      CellStyle leftHeaderStyle = getLeftHeaderStyle(resultWorkbook);
      CellStyle novelCenterStyle = getCenterStyle(novelSpeciesWorkbook);
      Sheet summarySheet = resultWorkbook.createSheet("Summary");
      Sheet overviewSheet = resultWorkbook.createSheet("Overview");
      Sheet ionModeSheet = resultWorkbook.createSheet("Polarity");
      Sheet adductsSheet = resultWorkbook.createSheet("Adducts");
      Sheet mzSheet = resultWorkbook.createSheet("mz");
      Sheet rtSheet = resultWorkbook.createSheet("Retention times");
      Sheet novelSheet = novelSpeciesWorkbook.createSheet("Novel Species");

      int summaryRowCount = 0;
      int detailsRowCount = 0;
      int novelRowCount = 0;
      
      String novelNrHeader = "No.";
      int novelNrColumn = 0;
      String novelSpeciesHeader = "Lipid species";
      int novelSpeciesColumn = 1;
      String novelMolSpeciesHeader = "Lipid molecular species";
      int novelMolSpeciesColumn = 2;

      Row summaryRow = summarySheet.createRow(summaryRowCount);
      summaryRowCount++;
      Row overviewRow = overviewSheet.createRow(detailsRowCount);
      Row ionModeRow = ionModeSheet.createRow(detailsRowCount);
      Row adductsRow = adductsSheet.createRow(detailsRowCount);
      Row mzRow = mzSheet.createRow(detailsRowCount);
      Row rtRow = rtSheet.createRow(detailsRowCount);
      detailsRowCount++;
      Row novelRow = novelSheet.createRow(novelRowCount);
      novelRowCount++;
      novelRowCount++;
      
      String speciesHeader = "Species";
      int speciesColumn = 0;
      int longestSpecies = 0;
      String molSpeciesHeader = "Molecular Species";
      int molSpeciesColumn = 1;
      int longestMolSpecies = 0;
      int summaryInfoColumn = 0;
      
      String longestPlatform = "";
      int longestValueSize = 0;
      
      //write header information
      writeHeaderCell(speciesHeader,speciesColumn,overviewRow,headerStyle);
      writeHeaderCell(speciesHeader,speciesColumn,ionModeRow,headerStyle);
      writeHeaderCell(speciesHeader,speciesColumn,adductsRow,headerStyle);
      writeHeaderCell(speciesHeader,speciesColumn,mzRow,headerStyle);
      writeHeaderCell(speciesHeader,speciesColumn,rtRow,headerStyle);
      writeHeaderCell(molSpeciesHeader,molSpeciesColumn,overviewRow,headerStyle);
      writeHeaderCell(molSpeciesHeader,molSpeciesColumn,ionModeRow,headerStyle);
      writeHeaderCell(molSpeciesHeader,molSpeciesColumn,adductsRow,headerStyle);
      writeHeaderCell(molSpeciesHeader,molSpeciesColumn,mzRow,headerStyle);
      writeHeaderCell(molSpeciesHeader,molSpeciesColumn,rtRow,headerStyle);

      createCell(novelRow, novelCenterStyle, novelNrColumn, novelNrHeader);
      createCell(novelRow, novelCenterStyle, novelSpeciesColumn, novelSpeciesHeader);
      createCell(novelRow, novelCenterStyle, novelMolSpeciesColumn, novelMolSpeciesHeader);
      
      int count = 0;
      Hashtable<Integer,Integer> totalLipidSpeciesDetected = new Hashtable<Integer,Integer>();
      Hashtable<Integer,String> platformNames = new Hashtable<Integer,String>();
      Hashtable<Integer,Integer> totalLipidMolecularSpeciesDetected = new Hashtable<Integer,Integer>();
      for (String platform : resultDetails.keySet()){
        count++;
        if (platform.length()>longestPlatform.length()) longestPlatform = platform;
        writeHeaderCell(platform,summaryInfoColumn+count,summaryRow,headerStyle);
        writeHeaderCell(platform,molSpeciesColumn+count,overviewRow,headerStyle);
        writeHeaderCell(platform,molSpeciesColumn+count,ionModeRow,headerStyle);
        writeHeaderCell(platform,molSpeciesColumn+count,adductsRow,headerStyle);
        writeHeaderCell(platform,molSpeciesColumn+count,mzRow,headerStyle);
        writeHeaderCell(platform,molSpeciesColumn+count,rtRow,headerStyle);
        createCell(novelRow, novelCenterStyle, novelMolSpeciesColumn+count, platform);
        totalLipidSpeciesDetected.put(count, 0);
        platformNames.put(count, platform);
        totalLipidMolecularSpeciesDetected.put(count, 0);
      }
      
      LinkedHashMap<String,LinkedHashMap<String,String>> novelSpecies = FoundBiologicalSpecies.getNovelSpecies();
      
      Hashtable<String,Hashtable<String,Hashtable<Integer,Integer>>> novelSpeciesDetections = new Hashtable<String,Hashtable<String,Hashtable<Integer,Integer>>>();
      
      for (String lipidClass : classSequence){
        for (String species : analyteSequence.get(lipidClass)){
          //TODO: write the Excel file - ATTENTION: grouping of molecular species might be problematic,
          //since "lipid species" are there too - check if they have been found in which file!!!
          String speciesId = lipidClass+" "+species;
          count = 0;
          Hashtable<Integer,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>> hitsOfOneSpecies = new Hashtable<Integer,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>>();
          for (String platform : resultDetails.keySet()){
            count++;
            Hashtable<String,LinkedHashMap<String,Vector<SingleAdductIdentificationVO>>> results = resultDetails.get(platform);
            if (results.containsKey(speciesId)) hitsOfOneSpecies.put(count, results.get(speciesId));
          }
          if (hitsOfOneSpecies.size()==0)continue;

          //write a lipid species column
          overviewRow = overviewSheet.createRow(detailsRowCount);
          ionModeRow = ionModeSheet.createRow(detailsRowCount);
          adductsRow = adductsSheet.createRow(detailsRowCount);
          mzRow = mzSheet.createRow(detailsRowCount);
          rtRow = rtSheet.createRow(detailsRowCount);
          detailsRowCount++;
          this.createCell(overviewRow, null, speciesColumn, speciesId);
          if (speciesId.length()>longestSpecies) longestSpecies = speciesId.length();
          this.createCell(ionModeRow, null, speciesColumn, speciesId);
          this.createCell(adductsRow, null, speciesColumn, speciesId);
          this.createCell(mzRow, null, speciesColumn, speciesId);
          this.createCell(rtRow, null, speciesColumn, speciesId);
          

          // first key: the molecular species name; second key: the platform it is present
          Hashtable<String,Hashtable<Integer,Vector<SingleAdductIdentificationVO>>> molSpeciesGrouped = new Hashtable<String,Hashtable<Integer,Vector<SingleAdductIdentificationVO>>>();
          
          for (int i=1; i!=(resultDetails.size()+1);i++){
            if (!hitsOfOneSpecies.containsKey(i)){
              //this.createCell(overviewRow, centerStyle, molSpeciesColumn+i, "not detectable");
              continue;
            }
            
            //for novel species
            if(speciesId.equalsIgnoreCase("Cer 24:2")){
              Hashtable<String,Hashtable<Integer,Integer>> cerHits = new Hashtable<String,Hashtable<Integer,Integer>>();
              if (novelSpeciesDetections.containsKey("Cer")) cerHits = novelSpeciesDetections.get("Cer");
              Hashtable<Integer,Integer> detections = new Hashtable<Integer,Integer>();
              if (cerHits.containsKey("Cer 24:2")) detections = cerHits.get("Cer 24:2");
              detections.put(i, i);
              cerHits.put("Cer 24:2", detections);
              novelSpeciesDetections.put("Cer",cerHits);
            }

            
            Hashtable<Integer,Integer> positiveSpeciesIdentification = new Hashtable<Integer,Integer>();
            Hashtable<Integer,Integer> negativeSpeciesIdentification = new Hashtable<Integer,Integer>();
            LinkedHashMap<String,String> adductsSpecies = new LinkedHashMap<String,String>();
            Vector<Double> firstRts = new Vector<Double>();
            Vector<Double> secondRts = new Vector<Double>();
            Vector<Double> thirdRts = new Vector<Double>();

            totalLipidSpeciesDetected.put(i, (totalLipidSpeciesDetected.get(i)+1));
            LinkedHashMap<String,Vector<SingleAdductIdentificationVO>> molecularSpeciesDetails = hitsOfOneSpecies.get(i);
            
            for (String molecularSpecies : molecularSpeciesDetails.keySet()){
              Vector<SingleAdductIdentificationVO> oneMolSpecies = molecularSpeciesDetails.get(molecularSpecies);
              Vector<SingleAdductIdentificationVO> pureStructures = new Vector<SingleAdductIdentificationVO>();
              String molSpeciesId = null;
              for (SingleAdductIdentificationVO oneAdduct : oneMolSpecies){
                adductsSpecies.put(oneAdduct.adduct, oneAdduct.measuredMz);
                for (Integer oneKey : oneAdduct.annotatedStructures.keySet()){
                  if (oneAdduct.ionMode.equalsIgnoreCase("+")) positiveSpeciesIdentification.put(oneKey, oneKey);
                  else if (oneAdduct.ionMode.equalsIgnoreCase("-")) negativeSpeciesIdentification.put(oneKey, oneKey);
                  else throw new Exception("There is a case which does not belong to +/- ion mode: "+molecularSpecies+"_"+oneAdduct.adduct+" "+platformNames.get(i));
                }
                //for the various retention times
                firstRts.addAll(new Vector<Double>(oneAdduct.firstRetentionTimes.values()));
                secondRts.addAll(new Vector<Double>(oneAdduct.secondRetentionTimes.values()));
                thirdRts.addAll(new Vector<Double>(oneAdduct.thirdRetentionTimes.values()));
                
                //removes lipid species only identifications
                if (molecularSpecies.equalsIgnoreCase(speciesId))continue;
                SingleAdductIdentificationVO structures = new SingleAdductIdentificationVO(oneAdduct);
                for (Integer oneKey : oneAdduct.annotatedStructures.keySet()){
                  if (oneAdduct.annotatedStructures.get(oneKey).equalsIgnoreCase(speciesId)) continue;
                  structures.annotatedStructures.put(oneKey, oneAdduct.annotatedStructures.get(oneKey));
                  if (oneAdduct.firstRetentionTimes.containsKey(oneKey))
                    structures.firstRetentionTimes.put(oneKey, oneAdduct.firstRetentionTimes.get(oneKey));
                  if (oneAdduct.secondRetentionTimes.containsKey(oneKey))
                    structures.secondRetentionTimes.put(oneKey, oneAdduct.secondRetentionTimes.get(oneKey));
                  if (oneAdduct.thirdRetentionTimes.containsKey(oneKey))
                    structures.thirdRetentionTimes.put(oneKey, oneAdduct.thirdRetentionTimes.get(oneKey));
                }
                if (structures.molFrequency!=structures.annotatedStructures.size())
                  System.out.println("The amount of structures does not correspond to the mol frequency: "+platformNames.get(i)+" "+structures.molecularSpecies);
                String anId = structures.molecularSpecies;
//                if (lipidClass.equalsIgnoreCase("DG")){
//                  anId = anId.replaceAll("/", "_");
//                  if (anId.indexOf("_")==anId.lastIndexOf("_")) anId += "_-";
//                }
                molSpeciesId = anId;
                pureStructures.add(structures);
              }
              
              if (pureStructures.size()==0) continue;
              Hashtable<Integer,Vector<SingleAdductIdentificationVO>> belongToOneMolSpecies = new Hashtable<Integer,Vector<SingleAdductIdentificationVO>>();
              
              String woPos = molSpeciesId.replaceAll("/", "_");
              for (String oneMol : molSpeciesGrouped.keySet()){
                String otherWoPos = oneMol.replaceAll("/", "_");
                if (!StaticUtils.isAPermutedVersion(woPos.substring(lipidClass.length()+1),otherWoPos.substring(lipidClass.length()+1),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) continue;
                belongToOneMolSpecies = molSpeciesGrouped.get(oneMol);
                molSpeciesId = oneMol;
              }
              //System.out.println(molSpeciesId+"_"+i);
              if (belongToOneMolSpecies.containsKey(i)) throw new Exception ("There is somehting wrong with the grouping of "+platformNames.get(i)+" "+pureStructures.get(0).molecularSpecies+" for molecular species");
              belongToOneMolSpecies.put(i, pureStructures);
              molSpeciesGrouped.put(molSpeciesId, belongToOneMolSpecies);
              //TODO: write the Excel file - ATTENTION: grouping of molecular species might be problematic,
              //since "lipid species" are there too - check if they have been found in which file!!!
            }
            
            String frequencyString = "found in "+(positiveSpeciesIdentification.size()+negativeSpeciesIdentification.size())+" MS runs";
            this.createCell(overviewRow, centerStyle, molSpeciesColumn+i, frequencyString);
            if (frequencyString.length()>longestValueSize) longestValueSize = frequencyString.length();
            String polarity = "";
            if (positiveSpeciesIdentification.size()>0) polarity = "+";
            if (polarity.length()>0 && negativeSpeciesIdentification.size()>0) polarity += "/";
            if (negativeSpeciesIdentification.size()>0) polarity += "-";
            this.createCell(ionModeRow, centerStyle, molSpeciesColumn+i, polarity);
            String adducts = "";
            String mz = "";
            for (String adduct : adductsSpecies.keySet()){
              if (adducts.length()>0){
                adducts+="; ";
                mz += "; ";
              }
              adducts += adduct;
              mz += adductsSpecies.get(adduct);
            }
            this.createCell(adductsRow, centerStyle, molSpeciesColumn+i, adducts);
            if (adducts.length()>longestValueSize) longestValueSize = adducts.length();
            this.createCell(mzRow, centerStyle, molSpeciesColumn+i, mz);
            //if (mz.length()>longestValueSize) longestValueSize = mz.length();
            String rtString = "";
            if (firstRts.size()>0){
              rtString += Calculator.FormatNumberToString(Calculator.mean(firstRts),2d);
              if (firstRts.size()>1 && Calculator.stddeviation(firstRts)>0.8d)
                System.out.println("Pay attention! The species "+speciesId+" at "+platformNames.get(i)+" shows a high retention time deviation "+Calculator.stddeviation(firstRts));
            }
            if (secondRts.size()>0){
              if (rtString.length()>0) rtString += "/";
              rtString += Calculator.FormatNumberToString(Calculator.mean(secondRts),2d);
              if (secondRts.size()>1 && Calculator.stddeviation(secondRts)>0.8d)
                System.out.println("Pay attention! The species "+speciesId+" at "+platformNames.get(i)+" shows a high retention time deviation "+Calculator.stddeviation(secondRts));
            }
            if (thirdRts.size()>0){
              if (rtString.length()>0) rtString += "/";
              rtString += Calculator.FormatNumberToString(Calculator.mean(thirdRts),2d);
              if (thirdRts.size()>1 && Calculator.stddeviation(thirdRts)>0.8d)
                System.out.println("Pay attention! The species "+speciesId+" at "+platformNames.get(i)+" shows a high retention time deviation "+Calculator.stddeviation(thirdRts));
            }
            this.createCell(rtRow, centerStyle, molSpeciesColumn+i, rtString);            
          }
          
          //first sort the lipid molecular species by their number of detections
          Hashtable<String,Integer> numberOfDetections = new Hashtable<String,Integer>();
          Hashtable<String,Hashtable<Integer,Vector<SingleAdductIdentificationVO>>> molSpeciesCorrectId = new Hashtable<String,Hashtable<Integer,Vector<SingleAdductIdentificationVO>>>();
          int highestIdentification = 0;
          
          for (String aMolSpecies : molSpeciesGrouped.keySet()){
            Hashtable<String,Integer> assignedSns = new Hashtable<String,Integer>();
            //unAssignedSns must not contain more than one value!!!
            Hashtable<String,Integer> unAssignedSns = new Hashtable<String,Integer>();
            int identified = 0;
            Hashtable<Integer,Vector<SingleAdductIdentificationVO>> oneMolSpecies = molSpeciesGrouped.get(aMolSpecies);
            for (Integer key : oneMolSpecies.keySet()){
              Vector<SingleAdductIdentificationVO> adducts = oneMolSpecies.get(key);
              for (SingleAdductIdentificationVO adduct : adducts){
                for (String struct : adduct.annotatedStructures.values()){
                  String structure = struct;
                  if (lipidClass.equalsIgnoreCase("DG")){

                    if ((structure.indexOf("_")!=-1 && structure.indexOf("_")==structure.lastIndexOf("_"))||
                        (structure.indexOf("/")!=-1 && structure.indexOf("/")==structure.lastIndexOf("/"))){
                      structure = structure.replaceAll("/", "_");
                      structure += "_-";
                    }
                  }

                  identified++;
                  int foundSns = 0;
                  if (structure.indexOf("/")!=-1){
                    if (assignedSns.containsKey(structure)) foundSns = assignedSns.get(structure);
                    foundSns++;
                    assignedSns.put(structure, foundSns);
                  }else{
                    String otherStructure = structure;
                    if (unAssignedSns.size()>0) otherStructure = unAssignedSns.keySet().iterator().next();
                    if (!StaticUtils.isAPermutedVersion(structure.substring(lipidClass.length()+1), otherStructure.substring(lipidClass.length()+1),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
                      throw new Exception ("There cannot be more than one FA combination for an unassigned structure "+platformNames.get(key)+" "+structure+"for molecular species");
                    foundSns++;
                    unAssignedSns.put(otherStructure, foundSns);
                  }
                }
              }
            }
            boolean moreThanOneSn = true;
            int highestNumber = 0;
            String correctId = null;
            for (String sn : assignedSns.keySet()){
              if (assignedSns.get(sn)>highestNumber){
                moreThanOneSn = false;
                highestNumber = assignedSns.get(sn);
                correctId = sn;
              }else if (assignedSns.get(sn)==highestNumber){
                moreThanOneSn = true;
                correctId = null;
              }
            }
            if (correctId==null){
              String idToChange = null;
              if (unAssignedSns.size()>0) idToChange = unAssignedSns.keySet().iterator().next();
              else idToChange = assignedSns.keySet().iterator().next().replaceAll("/", "_");
              correctId = lipidClass+" "+StaticUtils.sortFASequenceUnassigned(idToChange.substring(lipidClass.length()+1),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
            }
            correctId = checkForKnownManualCorrections(correctId);
            numberOfDetections.put(correctId, identified);
            molSpeciesCorrectId.put(correctId, oneMolSpecies);
            if (identified>highestIdentification) highestIdentification =  identified;
          }
          
          //now write the individual molecular species
          for (int nrIdents=highestIdentification; nrIdents!=-1; nrIdents--){
            //this is for sorting
            for (String molId : molSpeciesCorrectId.keySet()){
              if (numberOfDetections.get(molId)!=nrIdents) continue;
              Hashtable<Integer,Vector<SingleAdductIdentificationVO>> oneMolSpeciesHits = molSpeciesCorrectId.get(molId);
              overviewRow = overviewSheet.createRow(detailsRowCount);
              ionModeRow = ionModeSheet.createRow(detailsRowCount);
              adductsRow = adductsSheet.createRow(detailsRowCount);
              mzRow = mzSheet.createRow(detailsRowCount);
              rtRow = rtSheet.createRow(detailsRowCount);
              detailsRowCount++;
              this.createCell(overviewRow, null, molSpeciesColumn, molId);
              this.createCell(ionModeRow, null, molSpeciesColumn, molId);
              this.createCell(adductsRow, null, molSpeciesColumn, molId);
              this.createCell(mzRow, null, molSpeciesColumn, molId);
              this.createCell(rtRow, null, molSpeciesColumn, molId);
              if (molId.length()>longestMolSpecies) longestMolSpecies = molId.length();
              
              for (int i=1; i!=(resultDetails.size()+1);i++){
                if (!oneMolSpeciesHits.containsKey(i)){
                  //this.createCell(overviewRow, centerStyle, molSpeciesColumn+i, "not detectable");
                  continue;
                }
                Hashtable<Integer,Integer> positiveSpeciesIdentification = new Hashtable<Integer,Integer>();
                Hashtable<Integer,Integer> negativeSpeciesIdentification = new Hashtable<Integer,Integer>();
                LinkedHashMap<String,String> adductsSpecies = new LinkedHashMap<String,String>();
                Vector<Double> firstRts = new Vector<Double>();
                Vector<Double> secondRts = new Vector<Double>();
                Vector<Double> thirdRts = new Vector<Double>();

                totalLipidMolecularSpeciesDetected.put(i, (totalLipidMolecularSpeciesDetected.get(i)+1));
                Vector<SingleAdductIdentificationVO> singleAdducts = oneMolSpeciesHits.get(i);
                
                for (SingleAdductIdentificationVO oneAdduct : singleAdducts){
                  adductsSpecies.put(oneAdduct.adduct, oneAdduct.measuredMz);
                  for (Integer oneKey : oneAdduct.annotatedStructures.keySet()){
                    if (oneAdduct.ionMode.equalsIgnoreCase("+")) positiveSpeciesIdentification.put(oneKey, oneKey);
                    else if (oneAdduct.ionMode.equalsIgnoreCase("-")) negativeSpeciesIdentification.put(oneKey, oneKey);
                    else throw new Exception("There is a case which does not belong to +/- ion mode: "+molId+"_"+oneAdduct.adduct+" "+platformNames.get(i));
                  }
                  //for the various retention times
                  firstRts.addAll(new Vector<Double>(oneAdduct.firstRetentionTimes.values()));
                  secondRts.addAll(new Vector<Double>(oneAdduct.secondRetentionTimes.values()));
                  thirdRts.addAll(new Vector<Double>(oneAdduct.thirdRetentionTimes.values()));

                }
                String frequencyString = "found in "+(positiveSpeciesIdentification.size()+negativeSpeciesIdentification.size())+" MS runs";
                this.createCell(overviewRow, centerStyle, molSpeciesColumn+i, frequencyString);
                String polarity = "";
                if (positiveSpeciesIdentification.size()>0) polarity = "+";
                if (polarity.length()>0 && negativeSpeciesIdentification.size()>0) polarity += "/";
                if (negativeSpeciesIdentification.size()>0) polarity += "-";
                this.createCell(ionModeRow, centerStyle, molSpeciesColumn+i, polarity);
                String adducts = "";
                String mz = "";
                for (String adduct : adductsSpecies.keySet()){
                  if (adducts.length()>0){
                    adducts+="; ";
                    mz += "; ";
                  }
                  adducts += adduct;
                  mz += adductsSpecies.get(adduct);
                }
                this.createCell(adductsRow, centerStyle, molSpeciesColumn+i, adducts);
                if (adducts.length()>longestValueSize) longestValueSize = adducts.length();
                this.createCell(mzRow, centerStyle, molSpeciesColumn+i, mz);
                String rtString = "";
                if (firstRts.size()>0){
                  rtString += Calculator.FormatNumberToString(Calculator.mean(firstRts),2d);
                  if (firstRts.size()>1 && Calculator.stddeviation(firstRts)>0.8d)
                    System.out.println("Pay attention! The species "+molId+" at "+platformNames.get(i)+" shows a high retention time deviation "+Calculator.stddeviation(firstRts));
                }
                if (secondRts.size()>0){
                  if (rtString.length()>0) rtString += "/";
                  rtString += Calculator.FormatNumberToString(Calculator.mean(secondRts),2d);
                  if (secondRts.size()>1 && Calculator.stddeviation(secondRts)>0.8d)
                    System.out.println("Pay attention! The species "+molId+" at "+platformNames.get(i)+" shows a high retention time deviation "+Calculator.stddeviation(secondRts));
                }
                if (thirdRts.size()>0){
                  if (rtString.length()>0) rtString += "/";
                  rtString += Calculator.FormatNumberToString(Calculator.mean(thirdRts),2d);
                  if (thirdRts.size()>1 && Calculator.stddeviation(thirdRts)>0.8d)
                    System.out.println("Pay attention! The species "+molId+" at "+platformNames.get(i)+" shows a high retention time deviation "+Calculator.stddeviation(thirdRts));
                }
                this.createCell(rtRow, centerStyle, molSpeciesColumn+i, rtString);            

                //for novel species
                if (!novelSpecies.containsKey(lipidClass)) continue;
                LinkedHashMap<String,String> novMols = novelSpecies.get(lipidClass);
                String two = molId.replaceAll("/", "_").substring(lipidClass.length()+1);
                String rightOne = null;
                for (String novMol : novMols.keySet()){
                  String one = novMol.replaceAll("_0:0", "_-").replaceAll("\\(", "").replaceAll("\\)", "").replaceAll("/", "_");
                  if (lipidClass.equalsIgnoreCase("P-PE")) one = one.substring("PE".length()+1);
                  else one = one.substring(lipidClass.length()+1);
                  if (StaticUtils.isAPermutedVersion(one, two,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
                    if (lipidClass.equalsIgnoreCase("DG")){
                      for (SingleAdductIdentificationVO oneAdduct : singleAdducts){
                        for (String structure : oneAdduct.annotatedStructures.values()){
                          if (structure.indexOf("/")==-1) continue;
                          if (!oneAdduct.adduct.equalsIgnoreCase("+NH4+") && novMol.substring(novMol.lastIndexOf("/")).equalsIgnoreCase(structure.substring(structure.lastIndexOf("/")))){
                            rightOne = novMol;
                            break;
                          }
                        }
                        if (rightOne!=null)break;
                      }
                    }else{
                      rightOne = novMol;
                      break;
                    }
                  }
                }
                if (rightOne==null) continue;
                Hashtable<String,Hashtable<Integer,Integer>> novelHits = new Hashtable<String,Hashtable<Integer,Integer>>();
                if (novelSpeciesDetections.containsKey(lipidClass)) novelHits = novelSpeciesDetections.get(lipidClass);
                Hashtable<Integer,Integer> detections = new Hashtable<Integer,Integer>();
                if (novelHits.containsKey(rightOne)) detections = novelHits.get(rightOne);
                detections.put(i, i);
                novelHits.put(rightOne, detections);
                novelSpeciesDetections.put(lipidClass,novelHits);
              }
            }
          }
        }
      }
      
      Row speciesRow = summarySheet.createRow(summaryRowCount);
      writeHeaderCell("Lipid species identified",summaryInfoColumn,speciesRow,leftHeaderStyle);
      summaryRowCount++;
      Row molecularRow = summarySheet.createRow(summaryRowCount);
      writeHeaderCell("Lipid molecular species identified",summaryInfoColumn,molecularRow,leftHeaderStyle);
      summaryRowCount++;
      for (int i=1; i!=(resultDetails.size()+1);i++){
        this.createCell(speciesRow, centerStyle, summaryInfoColumn+i, String.valueOf(totalLipidSpeciesDetected.get(i)));
        this.createCell(molecularRow, centerStyle, summaryInfoColumn+i, String.valueOf(totalLipidMolecularSpeciesDetected.get(i)));
      }
      
      setColumnWidth(summarySheet, speciesColumn, "Lipid species identified", 0);
      setColumnWidth(overviewSheet, speciesColumn, speciesHeader, longestSpecies);
      setColumnWidth(ionModeSheet, speciesColumn, speciesHeader, longestSpecies);
      setColumnWidth(adductsSheet, speciesColumn, speciesHeader, longestSpecies);
      setColumnWidth(mzSheet, speciesColumn, speciesHeader, longestSpecies);
      setColumnWidth(rtSheet, speciesColumn, speciesHeader, longestSpecies);
      setColumnWidth(novelSheet, novelSpeciesColumn, "", longestSpecies);
      setColumnWidth(overviewSheet, molSpeciesColumn, molSpeciesHeader, longestMolSpecies);
      setColumnWidth(ionModeSheet, molSpeciesColumn, molSpeciesHeader, longestMolSpecies);
      setColumnWidth(adductsSheet, molSpeciesColumn, molSpeciesHeader, longestMolSpecies);
      setColumnWidth(mzSheet, molSpeciesColumn, molSpeciesHeader, longestMolSpecies);
      setColumnWidth(rtSheet, molSpeciesColumn, molSpeciesHeader, longestMolSpecies);
      setColumnWidth(novelSheet, novelMolSpeciesColumn, "", longestMolSpecies);
      
      count = 0;
      for (String platform : resultDetails.keySet()){
        count++;
        setColumnWidth(summarySheet, summaryInfoColumn+count, longestPlatform, 0);
        setColumnWidth(overviewSheet, molSpeciesColumn+count, longestPlatform, longestValueSize);
        setColumnWidth(ionModeSheet, molSpeciesColumn+count, longestPlatform, 0);
        setColumnWidth(adductsSheet, molSpeciesColumn+count, longestPlatform, longestValueSize);
        setColumnWidth(mzSheet, molSpeciesColumn+count, longestPlatform, longestValueSize);
        setColumnWidth(rtSheet, molSpeciesColumn+count, longestPlatform, longestValueSize);
        setColumnWidth(novelSheet, novelMolSpeciesColumn+count, "", longestPlatform.length());
      }

      int novelNr = 0;
      
      for (String lipidClass : novelSpecies.keySet()){
        LinkedHashMap<String,String> speciesNames = novelSpecies.get(lipidClass);
        if (!novelSpeciesDetections.containsKey(lipidClass)) throw new Exception("Novel species: There is nothing found for "+lipidClass);
        Hashtable<String,Hashtable<Integer,Integer>> classDetections = novelSpeciesDetections.get(lipidClass);
        String lastSpecies = "";
        for (String molSpecies : speciesNames.keySet()){
          String species = speciesNames.get(molSpecies);
          if (!classDetections.containsKey(molSpecies)) throw new Exception("Novel species: There is nothing found for "+molSpecies);
          Hashtable<Integer,Integer> detections = classDetections.get(molSpecies);
          if (detections.size()==0) throw new Exception("Novel species: There is nothing found for "+molSpecies);
          novelNr++;
          novelRow = novelSheet.createRow(novelRowCount);
          novelRowCount++;
          this.createCell(novelRow, novelCenterStyle, novelNrColumn, String.valueOf(novelNr));
          if (!species.equalsIgnoreCase(lastSpecies))
            this.createCell(novelRow, novelCenterStyle, novelSpeciesColumn, species);
          if (!lipidClass.equalsIgnoreCase("Cer"))
            this.createCell(novelRow, novelCenterStyle, novelMolSpeciesColumn, molSpecies);
          for (int i=1; i!=(resultDetails.size()+1);i++){
            if (detections.containsKey(i))
              this.createCell(novelRow, novelCenterStyle, novelMolSpeciesColumn+i, "yes");
          }
          
          lastSpecies = species;
        }
        novelRowCount++;
      }
      
      out = new BufferedOutputStream(new FileOutputStream(fileName));
      resultWorkbook.write(out);
      resultWorkbook.close();
      outNovel = new BufferedOutputStream(new FileOutputStream(fileNameNovel));
      novelSpeciesWorkbook.write(outNovel);
      outNovel.close();
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } finally {
      try {
        if (out!=null)out.close();
      }
      catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      try {
        if (outNovel!=null)outNovel.close();
      }
      catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }

    }
  }
  
  private void writeHeaderCell(String value, int column, Row row, CellStyle headerStyle){
    Cell cell = row.createCell(column);
    cell.setCellStyle(headerStyle);
    cell.setCellValue(value);
  }
  
  private String checkForKnownManualCorrections(String id) throws LipidCombinameEncodingException{
    if (id.startsWith("PC")){
      String fas = id.substring("PC ".length()).replaceAll("/", "_");
      if (StaticUtils.isAPermutedVersion(fas, "16:1_18:3",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PC 16:1_18:3";
      else if (StaticUtils.isAPermutedVersion(fas, "15:0_20:3",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PC 15:0_20:3";
      else if (StaticUtils.isAPermutedVersion(fas, "16:0_20:1",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PC 16:0_20:1";
      else if (StaticUtils.isAPermutedVersion(fas, "16:1_20:5",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PC 16:1_20:5";
      else if (StaticUtils.isAPermutedVersion(fas, "18:0_20:6",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PC 18:0_20:6";
      else if (StaticUtils.isAPermutedVersion(fas, "17:1_22:6",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PC 17:1_22:6";
      else if (StaticUtils.isAPermutedVersion(fas, "18:1_22:5",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PC 18:1_22:5";
      else if (StaticUtils.isAPermutedVersion(fas, "20:4_20:5",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PC 20:4_20:5";     
    }else if (id.startsWith("PE")){
      String fas = id.substring("PE ".length()).replaceAll("/", "_");
      if (StaticUtils.isAPermutedVersion(fas, "18:1_21:6",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PE 18:1_21:6";
    }else if (id.startsWith("PG")){
      String fas = id.substring("PG ".length()).replaceAll("/", "_");
      if (StaticUtils.isAPermutedVersion(fas, "18:2_22:6",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PG 18:2_22:6";
      else if (StaticUtils.isAPermutedVersion(fas, "20:4_22:6",LipidomicsConstants.CHAIN_SEPARATOR_NO_POS))
        return "PG 20:4_22:6";
    }
    return id;
  }
  
  private LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> getPositiveListOfExp1Standards() {
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> results = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>>();
    
    addStandardsToPositiveList("PI",results,getPIStandards());
    addStandardsToPositiveList("P-PC",results,getPPCStandards());
    addStandardsToPositiveList("P-PE",results,getPPEStandards());
    addStandardsToPositiveList("LPC",results,getLPCStandards());
    addStandardsToPositiveList("LPE",results,getLPEStandards());
    addStandardsToPositiveList("PS",results,getPSStandards());
    addStandardsToPositiveList("LPS",results,getLPSStandards());
    addStandardsToPositiveList("PC",results,getPCStandards());
    addStandardsToPositiveList("PE",results,getPEStandards());
    addStandardsToPositiveList("PG",results,getPGStandards());
    addStandardsToPositiveList("DG",results,getDGStandards());
    addStandardsToPositiveList("TG",results,getTGStandards());
    addStandardsToPositiveList("SM",results,getSMStandards());
    addStandardsToPositiveList("Cer",results,getCerStandards());  
    return results;
  }

    
  private void addStandardsToPositiveList(String lipidClass, LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> results,
      LinkedHashMap<String,ReferenceInfoVO> standards){
    LinkedHashMap<String,LinkedHashMap<String,String>> positives = new LinkedHashMap<String,LinkedHashMap<String,String>>();
    for (String species : standards.keySet()){
      String speciesString = new String(species);
      if (lipidClass.equalsIgnoreCase("TG") && speciesString.startsWith("d")) speciesString = speciesString.substring(1);
      LinkedHashMap<String,String> mols = new LinkedHashMap<String,String>();
      if (!lipidClass.equalsIgnoreCase("LPC") && !lipidClass.equalsIgnoreCase("LPE") && !lipidClass.equalsIgnoreCase("LPS") &&
          !lipidClass.equalsIgnoreCase("SM") && !lipidClass.equalsIgnoreCase("Cer")){
        ReferenceInfoVO info = standards.get(species);
        mols.put(info.getMS2Name(), info.getMS2Name());
      }
      positives.put(speciesString, mols);
    }
    results.put(lipidClass, positives);
  }
    
  private void positionalIsomersSmallerRange(){
    
//    String chromFile = "D:\\positionIsomers\\preExperiment\\Mix_5uM.chrom";
//    String excelFile = "D:\\positionIsomers\\preExperiment\\Mix_5uM_DG_34_1.xlsx";

    
//    String chromFile = "D:\\positionIsomers\\preExperiment\\Mix_5uM-2.chrom";
//    String excelFile = "D:\\positionIsomers\\preExperiment\\Mix_5uM-2_DG_34_1.xlsx";
//
//    
    String chromFile = "D:\\positionIsomers\\preExperiment\\Mix_5uM-95min.chrom";
    String excelFile = "D:\\positionIsomers\\preExperiment\\Mix_5uM-95min_DG_34_1.xlsx";
    try{
      
      String[] chromPaths = StringUtils.getChromFilePaths(chromFile);
      LipidomicsAnalyzer lAnalyzer;
      ElementConfigParser aaParser = Settings.getElementParser();
      aaParser.parse();
      lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);

      Vector<LipidParameterSet> hits = LDAResultReader.readResultFile(excelFile, new Hashtable<String,Boolean>()).getIdentifications().get("DG");
      for (LipidParameterSet setOld : hits){
        
        for (CgProbe probe : setOld.getIsotopicProbes().get(0)){
//          probe.LowerValley = probe.Peak-10f;
//          probe.UpperValley = probe.Peak+107f;
          
//          if (57.16<probe.Peak/60f&&probe.Peak/60f<57.26){
//            probe.LowerValley = probe.Peak-5f;
//            probe.UpperValley = probe.Peak+5f;
//          }else if (57.51<probe.Peak/60f&&probe.Peak/60f<57.61){
//            System.out.println("2.");
//          }
          
          if (57.16f<probe.Peak/60f||(68.8f<probe.Peak/60f && probe.Peak/60f<68.9f)){
            probe.UpperValley = probe.Peak+2f;
          }else if ((57.5f<probe.Peak/60f && probe.Peak/60f<57.7f) || 69.2f<probe.Peak/60f){
            probe.LowerValley = probe.Peak-2f;
          }
          
        }
        LipidParameterSet set = new LipidParameterSet(setOld); 
        MSnAnalyzer msnAnalyzer = new MSnAnalyzer("DG","Na",set,lAnalyzer,null,true,false);  
        LipidomicsMSnSet msn = (LipidomicsMSnSet)msnAnalyzer.getResult();
        
        for (Object nameObject : msn.getMSnIdentificationNames()){
          if (nameObject instanceof Vector){
            for (String name : (Vector<String>)nameObject){
              System.out.println(msn.getRt()+": "+name);
            }
          }else{
            String name = (String)nameObject;
            System.out.println(msn.getRt()+": "+name);
          }
        }
      }
      
    }catch (Exception ex){
      ex.printStackTrace();
    }
    
  }

  private LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> getPositiveListOfExp2Standards() {
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> results = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>>();
    
    addStandardsToPositiveList("LPC",results,getLPCExp2Standards());
    addStandardsToPositiveList("LPE",results,getLPEExp2Standards());
    addStandardsToPositiveList("PC",results,getPCExp2Standards());
    addStandardsToPositiveList("PE",results,getPEExp2Standards());
    return results;
  }

  private LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> getPositiveListOfExp3Standards() {
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> results = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>>();

    LinkedHashMap<String,LinkedHashMap<String,String>> analytes = new LinkedHashMap<String,LinkedHashMap<String,String>>();
    LinkedHashMap<String,String> structures = new LinkedHashMap<String,String>();
    structures.put("14:0/18:0", "14:0/18:0");
    structures.put("16:0/16:0", "16:0/16:0");
    analytes.put("32:0", structures);
    structures = new LinkedHashMap<String,String>();
    structures.put("16:0/18:0", "16:0/18:0");
    structures.put("17:0/17:0", "17:0/17:0");
    analytes.put("34:0", structures);
    structures = new LinkedHashMap<String,String>();
    structures.put("18:0/18:2", "18:0/18:2");
    structures.put("18:1/18:1", "18:1/18:1");
    analytes.put("36:2", structures);
    results.put("PC", analytes);

    analytes = new LinkedHashMap<String,LinkedHashMap<String,String>>();
    structures = new LinkedHashMap<String,String>();
    structures.put("18:0/18:2", "18:0/18:2");
    structures.put("18:1/18:1", "18:1/18:1");
    analytes.put("36:2", structures);
    structures = new LinkedHashMap<String,String>();
    structures.put("16:0/20:4", "16:0/20:4");
    structures.put("18:2/18:2", "18:2/18:2");
    analytes.put("36:4", structures);
    results.put("PE", analytes);

    analytes = new LinkedHashMap<String,LinkedHashMap<String,String>>();
    structures = new LinkedHashMap<String,String>();
    structures.put("18:1/18:1", "18:1/18:1");
    structures.put("18:0/18:2", "18:0/18:2");
    analytes.put("36:2", structures);
    results.put("PG", analytes);

    analytes = new LinkedHashMap<String,LinkedHashMap<String,String>>();
    structures = new LinkedHashMap<String,String>();
    structures.put("18:0/18:2", "18:0/18:2");
    structures.put("18:1/18:1", "18:1/18:1");
    analytes.put("36:2", structures);
    results.put("PS", analytes);
    
    analytes = new LinkedHashMap<String,LinkedHashMap<String,String>>();
    structures = new LinkedHashMap<String,String>();
    structures.put("16:0/18:0/16:0", "16:0/18:0/16:0");
    structures.put("19:0/12:0/19:0", "19:0/12:0/19:0");
    analytes.put("50:0", structures);
    results.put("TG", analytes);
    
     return results;
  }

  private float getCorrespondingExperiment3Area(String fileName, String machineName, String className, String molName,
      String structure, LinkedHashMap<String,String> posStructures, LipidomicsMSnSet set){
    float area = -1f;
    if (className.equalsIgnoreCase("PE") && molName.equalsIgnoreCase("36:4")){
      float rt = Float.parseFloat(set.getRt());
      if (machineName.equalsIgnoreCase("Orbitrap Velos Pro CID")||machineName.equalsIgnoreCase("Orbitrap Velos Pro HCD")){
        if ((structure.equalsIgnoreCase("18:2/18:2") && 22.6f<rt && rt<24.3f)||
          (structure.equalsIgnoreCase("16:0/20:4") && 24.3f<rt && rt<25.8f))
          area = set.getArea();
      } else if (machineName.equalsIgnoreCase("4000 QTRAP")){
        if (fileName.equalsIgnoreCase("Data20150826_Ex3_QTrap_neg-003_QTrap_Ex3.1_neg_Ex3_neg.xlsx")||
            fileName.equalsIgnoreCase("Data20150826_Ex3_QTrap_neg-015_QTrap_Ex3.3_neg_Ex3_neg.xlsx")||
            fileName.equalsIgnoreCase("Data20150826_Ex3_QTrap_neg-025_QTrap_Ex3.5_neg_Ex3_neg.xlsx")||
            fileName.equalsIgnoreCase("Data20150826_Ex3_QTrap_neg-029_QTrap_Ex3.5_neg_Ex3_neg.xlsx")){
          String[] myFAs = getUniqueFas(structure);
          String otherStructure = "";
          for (String struct : posStructures.keySet()){
            if (!struct.equalsIgnoreCase(structure)){
              otherStructure = struct;
              break;
            }
          }
          String[] otherFAs =  getUniqueFas(otherStructure);
          float[] areas = separateMSnPeak(set, myFAs, otherFAs);
          if (areas[0]>0 && areas[1]>0){
            area = areas[0];
          }
        }else if ((structure.equalsIgnoreCase("18:2/18:2") && 23.6f<rt && rt<24.6f)||
            (structure.equalsIgnoreCase("16:0/20:4") && 24.6f<rt && rt<25.8f))
            area = set.getArea();
      } else if (machineName.equalsIgnoreCase("SYNAPT G1 HDMS QTOF")){
        if ((structure.equalsIgnoreCase("18:2/18:2") && 17.3f<rt && rt<18.3f)||
            (structure.equalsIgnoreCase("16:0/20:4") && 18.3f<rt && rt<19.3f))
            area = set.getArea();
      }  
    }else{
      String[] myFAs = getUniqueFas(structure);
      String otherStructure = "";
      for (String struct : posStructures.keySet()){
        if (!struct.equalsIgnoreCase(structure)){
          otherStructure = struct;
          break;
        }
      }
      String[] otherFAs =  getUniqueFas(otherStructure);
      float[] areas = separateMSnPeak(set, myFAs, otherFAs);
      // For Orbitrap CID and 4000 QTRAP, in a second step, this has to be changed to "areas[0]>0"  to detect the ones, where only one species is present
      // it is no recommended to do this in the first run - then the FP peaks are in
      if (areas[0]>0 && areas[1]>0){
        area = areas[0];
      } else if (areas[0]>0 && machineName.equalsIgnoreCase("Orbitrap Velos Pro HCD")){  
        area = areas[0];
      } else if (machineName.equalsIgnoreCase("SYNAPT G1 HDMS QTOF")){
        if (className.equalsIgnoreCase("PG")||
            ((className.equalsIgnoreCase("PC")||className.equalsIgnoreCase("PE"))&&molName.equalsIgnoreCase("36:2"))){
          if (areas[0]>0)
            area = areas[0];
        }
      }
    }
    return area;
  }
  
  private String[] getUniqueFas(String structure){
    String[] notUnique = structure.split("/");
    Hashtable<String,String> unique = new Hashtable<String,String>();
    for (String fa: notUnique){
      unique.put(fa, fa);
    }
    String[] uni = new String[unique.size()];
    int count = 0;
    for (String fa : unique.keySet()){
      uni[count] = fa;
      count++;
    }
    return uni;
  }
  
  private void compareDeviationFromTheoreticalValueBasedOnIntensity(){
    //String file = "D:\\BiologicalExperiment\\Orbitrap_CID\\negative\\002_liver2-1_Orbitrap_CID_neg_negative.xlsx";
    //String file = "D:\\BiologicalExperiment\\Singapore\\negative\\Sample2.1_neg_01_negative.xlsx";
    String file = "D:\\Christer\\20170606\\NP MS1 R450k\\pos\\NP_MS1_pos_01_positive_LPC_PC_PE_TG.xlsx";
    try{
      ElementConfigParser aaParser = Settings.getElementParser();
//      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> idealMasses = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)QuantificationThread.parseQuantExcelFile("D:\\Christer\\20170531\\quant\\positive_LPC_PC_PE_TG.xlsx",  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f,false).get(3);
      aaParser.parse();
      Hashtable<String,Vector<LipidParameterSet>> classes = LDAResultReader.readResultFile(file,new Hashtable<String,Boolean>()).getIdentifications();
      List<IsotopicRatioDeviationVO> ratios = new ArrayList<IsotopicRatioDeviationVO>();
      for (String className : classes.keySet()){
        Vector<LipidParameterSet> analytes = classes.get(className);
//        Hashtable<String,Hashtable<String,QuantVO>> idealOfClass = idealMasses.get(className);
        for (LipidParameterSet set : analytes){
          ////if (!(set instanceof LipidomicsMSnSet)) continue;
          if (set.getIsotopicProbes().size()>1 && set.getIsotopicProbes().get(0).size()==1){
            float refArea = set.getIsotopicProbes().get(0).get(0).Area;
            Vector<Double> probabs = StaticUtils.calculateChemicalFormulaIntensityDistribution(aaParser, set.getChemicalFormula(), 5, false);
            int count = 0;          
/****            for (int i=1; i!=2set.getIsotopicProbes().size(); i++){
              Vector<CgProbe> probes = set.getIsotopicProbes().get(i);
              if (probes.size()!=1) continue;
              float area = probes.get(0).Area;
              float theoreticalRatio = (float)(probabs.get(i)/probabs.get(0));
              float ratio = area/refArea;
              IsotopicRatioDeviationVO devVO = new IsotopicRatioDeviationVO(className+set.getNameStringWithoutRt(),i,area,theoreticalRatio,
                  ratio);
              ratios.add(devVO);
            }*/ 
            Vector<CgProbe> probes = set.getIsotopicProbes().get(0);
//            QuantVO quant = idealOfClass.get(set.getNameStringWithoutRt()).get(set.getModificationName());
            for (CgProbe probe: probes){
              IsotopicRatioDeviationVO devVO = new IsotopicRatioDeviationVO(className+set.getNameString(),0,probe.Area,set.Mz[0],
                  (double)probe.Mz);
              ratios.add(devVO);
            }
          }
        }
      }
      Collections.sort(ratios,new GeneralComparator("at.tugraz.genome.TestClass$IsotopicRatioDeviationVO", "getArea", "java.lang.Float"));
      List<IsotopicRatioDeviationVO> ratiosDesc = new ArrayList<IsotopicRatioDeviationVO>();
      for (int i=(ratios.size()-1); i!=-1; i--){
        ratiosDesc.add(ratios.get(i));
      }
      
      // now plot the values to a chart
      ////XYSeries series1 = new XYSeries("Deviation from theoretical distribution");
/*      XYSeries series1 = new XYSeries("Deviation from theoretical m/z [mDa]");
      for (IsotopicRatioDeviationVO ratio : ratiosDesc){
        ////float percentDeviation = ((ratio.ratio_-ratio.theoreticalRatio_)*100f)/ratio.theoreticalRatio_;
        double percentDeviation = (ratio.getRatio()-ratio.getTheoreticalRatio())*1000d;
        String areaString = String.valueOf(ratio.area_);
        areaString.replaceAll("E", "e");
        ////System.out.println(areaString+"\t"+percentDeviation+"%\t"+ratio.name_+" +"+ratio.isotope_);
        System.out.println(areaString+"\t"+percentDeviation+"\t"+ratio.name_+" +"+ratio.isotope_);
        //XYWithName xyPoint = new XYWithName(ratio.getName(),(double)ratio.area_, (double)percentDeviation);
        //series1.add(xyPoint);
        series1.add((double)ratio.area_, (double)percentDeviation);
      }
      LogarithmicAxis xAxis = new LogarithmicAxis("Area under Curve");
      XYSeriesCollection coll1 = new XYSeriesCollection(series1);
      ////JFreeChart chart1 = ChartFactory.createScatterPlot("Ratio deviation", "Area under Curve", "Percentual deviation", coll1, PlotOrientation.VERTICAL, true, true, false);
      JFreeChart chart1 = ChartFactory.createScatterPlot("m/z deviation", "Area under Curve", "Delta m/z", coll1, PlotOrientation.VERTICAL, true, true, false);
      XYPlot plot = chart1.getXYPlot();
      plot.setDomainAxis(xAxis);
//      XYItemRenderer renderer = plot.getRenderer();
//      renderer.setToolTipGenerator(new TTGenerator());
      
      BufferedOutputStream stream;
      //stream = new BufferedOutputStream(new FileOutputStream("H:\\My Documents\\201706_Glucose_Paper_Martin\\DeviationInRelationToIntensity.png"));
      stream = new BufferedOutputStream(new FileOutputStream("D:\\Christer\\20170606\\NP MS1 R450k\\pos\\MzInRelationToIntensity.png"));
      ImageIO.write(chart1.createBufferedImage(1000, 700), "PNG", stream);
      stream.close();
*/
    }catch(Exception ex){
      ex.printStackTrace();
    }
  }
  
//  public class TTGenerator implements XYToolTipGenerator {
//
//    public String generateToolTip(XYDataset arg0, int arg1, int arg2)
//    {
//      if (arg0 instanceof XYWithName){
//        return ((XYWithName)arg0).getName();
//      }
//      return null;
//    }
//    
//  }
  
  public class IsotopicRatioDeviationVO{
    
    private String name_;
    private int isotope_;
    private Float area_;
    private double theoreticalRatio_;
    private double ratio_;
    
    
    
    public IsotopicRatioDeviationVO(String name, int isotope, Float area,
        double theoreticalRatio, double ratio)
    {
      super();
      this.name_ = name;
      this.isotope_ = isotope;
      this.area_ = area;
      this.theoreticalRatio_ = theoreticalRatio;
      this.ratio_ = ratio;
    }
    
    public String getName()
    {
      return name_;
    }
    public void setName(String name)
    {
      this.name_ = name;
    }
    public int getIsotope()
    {
      return isotope_;
    }
    public void setIsotope(int isotope)
    {
      this.isotope_ = isotope;
    }
    public Float getArea()
    {
      return area_;
    }
    public void setArea(Float area)
    {
      this.area_ = area;
    }
    public double getTheoreticalRatio()
    {
      return theoreticalRatio_;
    }
    public void setTheoreticalRatio(double theoreticalRatio)
    {
      this.theoreticalRatio_ = theoreticalRatio;
    }
    public double getRatio()
    {
      return ratio_;
    }
    public void setRatio(float ratio)
    {
      this.ratio_ = ratio;
    }
    
    
  }
  
//  public class XYWithName extends XYDataItem {
//    
//    private static final long serialVersionUID = 1L;
//    private String name_;
//    
//    public XYWithName(String name, double x, double y)
//    {
//      super(x, y);
//      this.name_ = name;
//    }
//
//    public String getName()
//    {
//      return name_;
//    }
//    
//    
//    
//  }

  private void readAlex123File(){
    String folderName = "D:\\Christer\\20170531\\target_lists_alex";
    try{
      TargetlistDirParser dirParser = new TargetlistDirParser(folderName,true);
      dirParser.parse();
      LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>>> sortedEntries = dirParser.getResults();
      for (String lipidClass : sortedEntries.keySet()){
        LinkedHashMap<String,LinkedHashMap<String,TargetlistEntry>> ofClass = sortedEntries.get(lipidClass);
        for (String analyte : ofClass.keySet()){
          LinkedHashMap<String,TargetlistEntry> ofAnalyte = ofClass.get(analyte);
          for (String mod : ofAnalyte.keySet()){
            TargetlistEntry entry = ofAnalyte.get(mod);
            //System.out.println(entry.getSpecies()+":\t"+entry.getAnalyteName()+" --------- "+entry.getDbs());
            System.out.println(entry.getSpecies()+":\t"+entry.getModName()+" -------- "+entry.getModFormula()+"\t"+entry.getAnalyteFormula());
          }
        }
//        System.out.println(lipidClass+": "+sortedEntries.get(lipidClass).size());
      }

    } catch (AlexTargetlistParserException ale){
      ale.printStackTrace();
    }
  }
  
  private void startQuantitationWithAlex123File(){
    //String chromFile = "D:\\Christer\\20170606\\NP MS1 R450k\\pos\\NP_MS1_pos_01.chrom";
    //String chromFile = "D:\\Christer\\20170606\\NP MS1 R450k\\neg\\NP_MS1_neg_01.chrom";
    String chromFile = "D:\\Christer\\20170310\\test\\002_liver2-1_Orbitrap_CID_pos.chrom";
    //String chromFile = "D:\\Christer\\20170310\\test2\\002_liver2-1_Orbitrap_CID_neg.chrom";
    String quantFile = "D:\\Christer\\20170531\\target_lists_alex";
    //String resultFile = "D:\\Christer\\20170606\\NP MS1 R450k\\pos\\NP_MS1_pos_01_target_lists_alex.xlsx";
    //String resultFile = "D:\\Christer\\20170606\\NP MS1 R450k\\neg\\NP_MS1_neg_01_target_lists_alex.xlsx";
    String resultFile = "D:\\Christer\\20170310\\test\\002_liver2-1_Orbitrap_CID_pos_target_lists_alex.xlsx";
    //String resultFile = "D:\\Christer\\20170310\\test2\\002_liver2-1_Orbitrap_CID_neg_target_lists_alex.xlsx";
    float cutoff = 0.02f;
    boolean positiveIonMode = true;
    try{
      LipidomicsConstants.getInstance().setRelativeMS1BasePeakCutoff(String.valueOf(cutoff));
      QuantificationThread quantThread = new QuantificationThread(chromFile,quantFile,resultFile, 
          5f, 5f, 2, 1, true, cutoff, 0f, 7,positiveIonMode,false);
      quantThread.start();
      while (!quantThread.finished()){
        try {
          this.wait(1000);
        }
        catch (Exception ex) {
        }
      }
      System.out.println("Quant is finished!");
      if (quantThread.getErrorString()!=null){
        System.out.println("ERROR: "+quantThread.getErrorString());
      }
    }catch(Exception ex){
      ex.printStackTrace();
    }
  }
  
  private void readIndexFile(){
    //String indexFilePath = "D:\\positionIsomers\\experiment\\Mix_5uM-95min-1.chrom\\Mix_5uM-95min-1.idx2";
//    String indexFilePath = "D:\\positionIsomers\\test\\Mix_5uM-95min-1.chrom\\Mix_5uM-95min-1.idx2";
    String indexFilePath = "D:\\Kristaps\\20171129\\TG quant NIST\\MCC007_Lipid01_NIST1_20171124_negative.chrom\\MCC007_Lipid01_NIST1_20171124_negative.idx2";
    DataInputStream inStream = null;
    try{
      inStream = new DataInputStream(new FileInputStream(indexFilePath));
      int count = 0;
      while (inStream.available()>0){
        int currentLine = inStream.readInt();
        long bytesToSkip = inStream.readLong();
        System.out.println(currentLine+": "+bytesToSkip);
//          cachedLine.put(new Integer(count), new Integer(currentLine));
//          cachedIndex.put(new Integer(count), new Long(bytesToSkip));
        count++;
      }
    }catch(Exception iox){
      iox.printStackTrace();;
    }finally{
      if (inStream!=null)
        try {inStream.close();}catch (IOException iox){iox.printStackTrace();}
    }
  }
  
  private void readRttFile(){
    String rttFilePath = "D:\\shotgun\\20180206\\050_Shotgun1_IS_negative.chrom\\050_Shotgun1_IS_negative.rtt";
    DataInputStream inStream = null;
    try{
      inStream = new DataInputStream(new FileInputStream(rttFilePath));
      int count = 0;
      while (inStream.available()>0){
        int currentLine = inStream.readInt();
        float time = inStream.readFloat();
        System.out.println(currentLine+": "+time);
//          cachedLine.put(new Integer(count), new Integer(currentLine));
//          cachedIndex.put(new Integer(count), new Long(bytesToSkip));
        count++;
      }
    }catch(Exception iox){
      iox.printStackTrace();;
    }finally{
      if (inStream!=null)
        try {inStream.close();}catch (IOException iox){iox.printStackTrace();}
    }
  }
  
  private void detectMsnByAlex123(){
    String chromFile = "D:\\Christer\\20170310\\test\\002_liver2-1_Orbitrap_CID_pos.chrom";
    String quantFile = "D:\\Christer\\20170531\\target_lists_alex";
    boolean positiveIonMode = true;
    try {
      Vector quants = QuantificationThread.getCorrectAnalyteSequence(quantFile, positiveIonMode);
      System.out.println("---------------------");
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantObjects = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)quants.get(3);
      TargetlistEntry quant = (TargetlistEntry)quantObjects.get("TAG").get("40:0").get("+NH4+");
//      Hashtable<Integer,Hashtable<String,TargetlistEntry>> fragments = quant.getMsnFragments();
//      for (Integer msLevel : fragments.keySet()){
//        System.out.println("msLevel: "+msLevel.toString());
//        Hashtable<String,TargetlistEntry> fragmentsOfLevel = fragments.get(msLevel);
//        for (String name : fragmentsOfLevel.keySet()){
//          System.out.println(name);
//        }
//      }
      String[] chromPaths = StringUtils.getChromFilePaths(chromFile);
      LipidomicsAnalyzer lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
      setStandardParameters(lAnalyzer);
      Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> isotopicProbes = lAnalyzer.processByMzProbabsAndPossibleRetentionTime((float)quant.getAnalyteMass(),quant.getCharge(),  
          -1f, -1f, -1f, -1, quant.getMustMatchProbabs(), quant.getProbabs(),1,false);
      for (Hashtable<Integer,Vector<CgProbe>> oneHit : isotopicProbes.values()){
        LipidParameterSet param = null;
        ////LipidParameterSet param = SingleQuantThread.createLipidParameterSet(oneHit,quant.getNegativeStartValue(), (float)quant.getAnalyteMass(),
        ////    quant.getAnalyteName(), quant.getDbs(), quant.getModName(), quant.getAnalyteFormula(), quant.getModFormula(), 
        ////    quant.getCharge());
        System.out.println("RT: "+param.getRt()+";"+(float)quant.getAnalyteMass());
        //MSnAnalyzer analyzer = new MSnAnalyzer(null,"TG","NH4",param,lAnalyzer,quant,false,false);
        MSnAnalyzer analyzer = new MSnAnalyzer(null,quant.getAnalyteClass(),quant.getModName(),param,lAnalyzer,quant,false,true,false);
        if (!(analyzer.getResult() instanceof LipidomicsMSnSet)) continue;
        LipidomicsMSnSet result = (LipidomicsMSnSet)analyzer.getResult();
        for (Object nameObject : result.getMSnIdentificationNames()){
          if (nameObject instanceof Vector){
            for (String name : (Vector<String>)nameObject){
              System.out.println(result.getRt()+": "+name);
            }
          }else{
            String name = (String)nameObject;
            System.out.println(result.getRt()+": "+name);
          }
        }
      }
    }
    catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  private int readIndexFile(DataInputStream inStream, Hashtable<Integer,Integer> lines, Hashtable<Integer,Long> indices,
      int count, long previousLength) throws IOException{
    int iteration = 0;
    while (inStream.available()>0){
      int currentLine = inStream.readInt();
      long bytesToSkip = inStream.readLong();
      //fill empty spaces at end of each by thread translated file
      if (iteration==0 && count>0){
        while ((lines.get(count-1)+CgDefines.numberOfEntriesForIndex)<currentLine){
          lines.put(count, lines.get(count-1)+CgDefines.numberOfEntriesForIndex);
          indices.put(count, previousLength);
          count++;
        }
      }
      lines.put(new Integer(count), new Integer(currentLine));
      indices.put(new Integer(count), (new Long(bytesToSkip))+previousLength);
      count++;
      iteration++;
    }
    return count;
  }
  
  private void mergeIdx2(){
//    try{
//      //DataInputStream inStream = new DataInputStream(new FileInputStream("D:\\Alex\\20170913\\data\\20170904-HILIC-Neg-MSMS-03-TQC01_working.chrom\\20170904-HILIC-Neg-MSMS-03-TQC01.idx2"));
//      DataInputStream inStream = new DataInputStream(new FileInputStream("D:\\Alex\\20170913\\data\\20170904-HILIC-Neg-MSMS-03-TQC01.chrom\\13\\20170904-HILIC-Neg-MSMS-03-TQC01.idx2"));
//      while (inStream.available()>0){
//        int currentLine = inStream.readInt();
//        long bytesToSkip = inStream.readLong();
//        System.out.println(currentLine+";"+bytesToSkip);
//      }
//    }catch(Exception ex){
//      ex.printStackTrace();
//    }
    try{
      int count = 0;
      long previousLength = 0;    
      Hashtable<Integer,Integer> lines = new Hashtable<Integer,Integer>();
      Hashtable<Integer,Long> indices = new Hashtable<Integer,Long>();

      DataInputStream indexStream = new DataInputStream(new FileInputStream("D:\\Alex\\20170913\\data\\20170904-HILIC-Neg-MSMS-03-TQC01_merge.chrom\\20170904-HILIC-Neg-MSMS-03-TQC01.idx2"));
      count = readIndexFile(indexStream,lines,indices,count,previousLength);
      indexStream.close();
      previousLength += (new File ("D:\\Alex\\20170913\\data\\20170904-HILIC-Neg-MSMS-03-TQC01_merge.chrom\\20170904-HILIC-Neg-MSMS-03-TQC01.chrom2")).length();

      String baseDir = "D:\\Alex\\20170913\\data\\20170904-HILIC-Neg-MSMS-03-TQC01_merge.chrom\\";
      for (int j=1; j!=28; j++){
        String dir = baseDir+String.valueOf(j)+"/";
        String indexFileName = dir+"20170904-HILIC-Neg-MSMS-03-TQC01"+".idx2";
        DataInputStream in = new DataInputStream(new FileInputStream(indexFileName));
        count = readIndexFile(in,lines,indices,count,previousLength);
        in.close();
        previousLength += (new File(dir+"20170904-HILIC-Neg-MSMS-03-TQC01"+".chrom2").length());
      }
      for (int i=0;i!=count;i++){
        System.out.println(lines.get(i)+";"+indices.get(i));
      }
    }catch(Exception ex){
      ex.printStackTrace();
    }
  }

  @Override
  public float getLowerThreshold()
  {
    // TODO Auto-generated method stub
    return 0;
  }

  @Override
  public float getUpperThreshold()
  {
    // TODO Auto-generated method stub
    return 0;
  }

  private void faSummaryConverter(){
    LDAToFASummaryConverter converter = new LDAToFASummaryConverter("D:\\Experiment1\\Orbitrap_CID\\positive\\50");
    try {
      converter.convert();
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  private void mzTabValidation(){
//    String filePath = "D:\\01-mztab.txt";
//    MZTabFileParser parser = new MZTabFileParser(new File(filePath));
//    try {
//      parser.parse(System.out);
//      MzTab mzTabFile = parser.getMZTabFile();
//      System.out.println("Software: "+mzTabFile.getMetadata().getSoftware().get(0).getParameter().getName());
//      System.out.println("11111111111111111");
//      MZTabFileValidator validator = new MZTabFileValidator(parser.getErrorList());
//      validator.toString();
//
//    }
//    catch (IOException e) {
//      // TODO Auto-generated catch block
//      e.printStackTrace();
//    }
    String filePath = "D:\\01-mztab.txt";
    //String filePath = "D:\\01-mztab_test.txt";
    //String filePath = "D:\\lda2-lipidomics.mztab";
    File inFile = new File(filePath);
    System.out.println(
        "Beginning validation of mztab file: " + inFile.
            getAbsolutePath());
    try {
      OutputStream out = System.out;
      MZTabErrorType.Level level = MZTabErrorType.Level.Error;
      MzTabFileParser mzTabParser = new MzTabFileParser(inFile);
      MZTabErrorList errorList = mzTabParser.parse(out, level);
      System.out.println("111111111111111111");
      if (!errorList.isEmpty()) {
//        System.out.println("There were errors while processing your file, please check the output for details!");
        List<ValidationMessage> messages = errorList.convertToValidationMessages();
        for (ValidationMessage message : messages){
          System.out.println("Message: "+message.getMessage());
        }
      }
    } catch (IOException e) {
      System.out.println(
          "Caught an IO Exception: " + e.getMessage());
    }
  }
  
  private void evaluatePrmData(){
    String chromFile = "C:\\data\\Alex\\PRM\\PRM_Beispiel__006_Platelets_PRM_CerOH_HexCerOH.chrom";
//    String quantFile = "D:\\Christer\\20170531\\target_lists_alex";
//    boolean positiveIonMode = true;
    float mz = 716.522481811612f;
    int charge = 1;
    int msLevel = 2;
    String className = "PE";
    String modName = "H";
    String analyteName = "34:2";
    String formula = "H75 P1 O8 N1 C39";
    try {
      String[] chromPaths = StringUtils.getChromFilePaths(chromFile);

      LipidomicsAnalyzer lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
      setStandardParameters(lAnalyzer);
      Hashtable<Integer,Hashtable<Integer,Vector<CgProbe>>> results = lAnalyzer.processPrmData(mz, charge, msLevel, className, modName, analyteName, formula,0);
      for (Integer hitNumber : results.keySet()){
        for (Integer isoNr: results.get(hitNumber).keySet()){
          for (CgProbe probe : results.get(hitNumber).get(isoNr)){
            System.out.println("! "+probe.Peak/60f+" "+probe.Area);
            System.out.println(probe.LowerValley/60f+"-"+probe.UpperValley/60f);
          }
        }
      }
      
      
    } catch (Exception ex){
      ex.printStackTrace();
    }
  }
  private void batchQuantByCommandLine() {
    //TODO: this should be input parameters
    String rawDirString = "/home/juergen/data/Orbitrap-CID/-50";
    String quantDirString = "/home/juergen/data/Orbitrap-CID/-50/massList";
    int amountOfIsotopes = 2;
    int isotopesMustMatch = 1;
    float minusTimeTol = 0f;
    float plusTimeTol = 0f;
    boolean searchUnknownBatchTime = true;
    float cutoff = 0.0001f;
    float rtShift = 0f;
    int nrProcessors = 11;
    
    BatchQuantificationTableModel batchQuantTableModel = new BatchQuantificationTableModel();
    BatchQuantificationTable batchQuantTable = new BatchQuantificationTable(batchQuantTableModel);
    
    //TODO: these are visual components that are not required for command line
    JLabel quantifyingBatchLabel = new JLabel("Quantifying");
    JProgressBar progressBatchBar = new JProgressBar();
    progressBatchBar.setMaximum(100);

    
    Vector<File> rawFiles = new Vector<File>();
    Vector<File> quantFiles = new Vector<File>();
    if (rawDirString!=null && rawDirString.length()>0 && quantDirString!=null&&quantDirString.length()>0){
      File rawDir = new File(rawDirString );
      File quantDir = new File(quantDirString);
      if (rawDir.exists()&&rawDir.isDirectory()&&quantDir.exists()&&quantDir.isDirectory()){
        File[] rawFileCandidates = rawDir.listFiles();
        Hashtable<String,Vector<File>> avoidDuplication = new Hashtable<String,Vector<File>>();
        boolean mzXMLOrChromPresent = false;
        for (int i=0; i!=rawFileCandidates.length;i++){
          if (rawFileCandidates[i].isFile()){
            String[] fileNameAndSuffix = StaticUtils.extractFileNameAndSuffix(rawFileCandidates[i].getAbsolutePath()); 
            String suffix = fileNameAndSuffix[1];
            String fileName = fileNameAndSuffix[0];
            if (suffix.equalsIgnoreCase("mzxml")||suffix.equalsIgnoreCase("raw")||suffix.equalsIgnoreCase("chrom")||suffix.equalsIgnoreCase("wiff")){
              if (suffix.equalsIgnoreCase("mzxml")||suffix.equalsIgnoreCase("chrom")) mzXMLOrChromPresent = true;
              Vector<File> theFiles = new Vector<File>();
              if (avoidDuplication.containsKey(fileName)){
                
                theFiles = avoidDuplication.get(fileName);
                }
                theFiles.add(rawFileCandidates[i]);
                avoidDuplication.put(fileName, theFiles);
              }
            }
            if (rawFileCandidates[i].isDirectory()){
              String[] fileNameAndSuffix = StaticUtils.extractFileNameAndSuffix(rawFileCandidates[i].getAbsolutePath()); 
              String suffix = fileNameAndSuffix[1];
              String fileName = fileNameAndSuffix[0];
              if (suffix.equalsIgnoreCase("raw")|| suffix.equalsIgnoreCase("d") ||suffix.equalsIgnoreCase("chrom")){
                if (suffix.equalsIgnoreCase("chrom")) mzXMLOrChromPresent = true;
                Vector<File> theFiles = new Vector<File>();
                if (avoidDuplication.containsKey(fileName)){
                  theFiles = avoidDuplication.get(fileName);
                }
                theFiles.add(rawFileCandidates[i]);
                avoidDuplication.put(fileName, theFiles);
              }
            }
          }
          for (String key : avoidDuplication.keySet()){
            Vector<File> theFiles = avoidDuplication.get(key);
            if (theFiles.size()==1){
              String suffix = StaticUtils.extractFileNameAndSuffix(theFiles.get(0).getAbsolutePath())[1];
              if (!mzXMLOrChromPresent || !suffix.equalsIgnoreCase("wiff"))
                rawFiles.add(theFiles.get(0));
            }else{
              int selectedIndex = -1;
              for (int i=0; i!=theFiles.size();i++){
                File file = theFiles.get(i);
                String suffix = file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf(".")+1);
                if (mzXMLOrChromPresent && suffix.equalsIgnoreCase("wiff")) continue;
                if (suffix.equalsIgnoreCase("chrom")){
                  selectedIndex = i;
                }
              }
              if (selectedIndex>-1){
                rawFiles.add(theFiles.get(selectedIndex));
              }else{
                for (int i=0; i!=theFiles.size();i++){
                  File file = theFiles.get(i);
                  String suffix = file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf(".")+1);
                  if (mzXMLOrChromPresent && suffix.equalsIgnoreCase("wiff")) continue;
                  if (suffix.equalsIgnoreCase("mzXML")){
                    rawFiles.add(theFiles.get(i));
                  }
                }  
              }
            }
          }
          File[] quantificationFileCandidates = quantDir.listFiles();
          boolean containsTxtFiles = false;
          for (int i=0; i!=quantificationFileCandidates.length;i++){
            String suffix = quantificationFileCandidates[i].getAbsolutePath().substring(quantificationFileCandidates[i].getAbsolutePath().lastIndexOf(".")+1);
            if (suffix.equalsIgnoreCase("xls")||suffix.equalsIgnoreCase("xlsx")){
              quantFiles.add(quantificationFileCandidates[i]);
            } else if (suffix.equalsIgnoreCase("txt"))
              containsTxtFiles = true;
          }
          if (Settings.useAlex() && containsTxtFiles)
            quantFiles.add(quantDir);
          Vector<RawQuantificationPairVO> pairs = new Vector<RawQuantificationPairVO>();
          
          if (rawFiles.size()>0 && quantFiles.size()>0 && (pairs = LipidDataAnalyzer.generateQuantificationPairVOs(rawFiles,quantFiles)).size()>0){
            boolean ionMode = false;
//            if (this.ionModeBatch_!=null && ((String)ionModeBatch_.getSelectedItem()).equalsIgnoreCase("+"))
//              ionMode = true;
            batchQuantTableModel.clearFiles();
            batchQuantTableModel.addFiles(pairs);
            progressBatchBar.setValue(0);
            BatchQuantThread batchQuantThread_ = new BatchQuantThread(batchQuantTable, batchQuantTableModel,progressBatchBar, 
              quantifyingBatchLabel, minusTimeTol,plusTimeTol,amountOfIsotopes,isotopesMustMatch,searchUnknownBatchTime, cutoff, 
                rtShift, nrProcessors,ionMode,true);
            batchQuantThread_.start();
          }else{
            if (rawFiles.size()==0){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "In the specified raw directory are no quantifyable files");
            }
            if (quantFiles.size()==0){
              @SuppressWarnings("unused")
              WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "In the specified quant directory are no quantifyable files");
            }
            if (rawFiles.size()>0 && quantFiles.size()>0)
              new WarningMessage(new JFrame(), "Warning", "In the specified directories are no quantifyable raw/quant pairs");
          }
        }else{
          if (!rawDir.exists()||!rawDir.isDirectory()){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The raw directory does not exist");
          }
          if (!quantDir.exists()||!quantDir.isDirectory()){
            @SuppressWarnings("unused")
            WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The quantification directory does not exist");
          }
        }

    }
  }
  
  private void groupAlexResultsByRt() {
    AlexRtGrouper grouper = new AlexRtGrouper("C:\\data\\Christer\\20190111\\positive",0.4d);
    try {
      grouper.groupTheEntries();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
  
  private void generateSphingolipidMassList() {
    try{
      ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
      parser.parse();
      
      //for sphingoid backbone chain library
      //String quantFile = "C:\\Development\\LipidDataAnalyzer\\fattyAcids\\dLCB.xlsx";
      
      //for sphingoid bases
      //String quantFile = "C:\\Sphingolipids\\massLists\\tSphBase_positive.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\dSphBase_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\tSphBase_positive.xlsx";
      
      //for Cer
      //String quantFile = "C:\\Sphingolipids\\massLists\\mCer_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\dCer_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\tCer_positive.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\qCer_positive.xlsx";
      
      //for Cer1P
      //String quantFile = "C:\\Sphingolipids\\massLists\\mCer1P_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\dCer1P_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\tCer1P_positive.xlsx";

      //for SM
      //String quantFile = "C:\\Sphingolipids\\massLists\\mSM_positive.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\dSM_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\tSM_negative.xlsx";

      //for LSM
      //String quantFile = "C:\\Sphingolipids\\massLists\\mLSM_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\dLSM_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\tLSM_negative.xlsx";

      //for HexCer
      //String quantFile = "C:\\Sphingolipids\\massLists\\mHexCer_negative.xlsx";
      String quantFile = "C:\\Sphingolipids\\massLists\\dHexCer_positive.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\tHexCer_negative.xlsx";

      //for LacCer
      //String quantFile = "C:\\Sphingolipids\\massLists\\dLacCer_negative.xlsx";
      
      //for S1P
      //String quantFile = "C:\\Sphingolipids\\massLists\\mS1P_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\dS1P_negative.xlsx";
      //String quantFile = "C:\\Sphingolipids\\massLists\\tS1P_negative.xlsx";
      
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(quantFile));
      XSSFWorkbook resultWorkbook = new XSSFWorkbook();
      XSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
      
      //for sphingoid backbone chain library
      //Sheet resultSheet = resultWorkbook.createSheet("dLCB");
      
      //for sphingoid bases
      //Sheet resultSheet = resultWorkbook.createSheet("mSphBase");
      //Sheet resultSheet = resultWorkbook.createSheet("dSphBase");
      //Sheet resultSheet = resultWorkbook.createSheet("tSphBase");

      //for Cer
      //Sheet resultSheet = resultWorkbook.createSheet("mCer");
      //Sheet resultSheet = resultWorkbook.createSheet("dCer");
      //Sheet resultSheet = resultWorkbook.createSheet("tCer");
      //Sheet resultSheet = resultWorkbook.createSheet("qCer");

      //for Cer1P
      //Sheet resultSheet = resultWorkbook.createSheet("mCer1P");
      //Sheet resultSheet = resultWorkbook.createSheet("dCer1P");
      //Sheet resultSheet = resultWorkbook.createSheet("tCer1P");

      //for SM
      //Sheet resultSheet = resultWorkbook.createSheet("mSM");
      //Sheet resultSheet = resultWorkbook.createSheet("dSM");
      //Sheet resultSheet = resultWorkbook.createSheet("tSM");

      //for LSM
      //Sheet resultSheet = resultWorkbook.createSheet("mLSM");
      //Sheet resultSheet = resultWorkbook.createSheet("dLSM");
      //Sheet resultSheet = resultWorkbook.createSheet("tLSM");

      //for HexCer
      //Sheet resultSheet = resultWorkbook.createSheet("mHexCer");
      Sheet resultSheet = resultWorkbook.createSheet("dHexCer");
      //Sheet resultSheet = resultWorkbook.createSheet("tHexCer");

      //for LacCer
      //Sheet resultSheet = resultWorkbook.createSheet("dLacCer");

      //for S1P
      //Sheet resultSheet = resultWorkbook.createSheet("mS1P");
      //Sheet resultSheet = resultWorkbook.createSheet("dS1P");
      //Sheet resultSheet = resultWorkbook.createSheet("tS1P");

      
      int rowCount = 4;
      Row outRow = resultSheet.createRow(rowCount);
      Cell label = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("Name");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("dbs");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("C");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("H");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("O");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("P");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("N");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(8,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("D");
      label.setCellStyle(headerStyle);
      //for sphingoid backbone chain library
//      label = outRow.createCell(9,HSSFCell.CELL_TYPE_STRING);
//      label.setCellValue("mass");
//      label.setCellStyle(headerStyle);
      //for all others
      label = outRow.createCell(9,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("M");
      label.setCellStyle(headerStyle);

      
      //for positive ion mode data
      label = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);      
      label.setCellValue("mass(form[+H] name[H])");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);      
      label.setCellValue("mass(form[-OH] name[-OH])");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(12,HSSFCell.CELL_TYPE_STRING);      
      label.setCellValue("mass(form[+Na] name[Na])");
      label.setCellStyle(headerStyle);
      //for negative ion mode data
//      label = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);      
//      label.setCellValue("mass(form[-H] name[-H])");
//      label.setCellStyle(headerStyle);
//      label = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);      
//      label.setCellValue("mass(form[+HCOO] name[HCOO])");
//      label.setCellStyle(headerStyle);
//      label = outRow.createCell(12,HSSFCell.CELL_TYPE_STRING);      
//      label.setCellValue("mass(form[+Cl] name[Cl])");
//      label.setCellStyle(headerStyle);
      
      //for all except for sphingoid backbone chain library
      label = outRow.createCell(13,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("tR (min)");
      label.setCellStyle(headerStyle);
      rowCount++;
      
      Hashtable<Integer,Integer> cDbsCombi = new Hashtable<Integer,Integer>();

      //for sphingoid backbone chain library
//      for (int i=10; i!=21; i++) {
//        for (int j=0; j!=5; j++) {
//          if (j>4 && i<12)
//            continue;
//          if (j>5 && i<17)
//            continue;
//          cDbsCombi.put(i, j);
//        }
//      }
      
      //for sphingoid bases, LSM and S1P
//      for (int i=2; i!=31; i++) {
//        for (int j=0; j!=6; j++) {
//          if (j>0 && i<10)
//            continue;
//          cDbsCombi.put(i, j);
//        }
//      }
      
      //for Cer, Cer1P, SM and HexCer
      for (int i=20; i!=49; i++) {
        for (int j=0; j!=8; j++) {
          if (j>1 && i<26)
            continue;
          else if (j>2 && i<28)
            continue;
          else if (j>3 && i<32)
            continue;
          else if (j>4 && i<34)
            continue;
          else if (j>5 && i<36)
            continue;
          cDbsCombi.put(i, j);
        }
      }
      
      //for sphingoid backbone chain library
      //for (int cAtoms=10; cAtoms<21; cAtoms+=1){
      
      //for sphingoid bases, LSM and S1P
      //for (int cAtoms=2; cAtoms<31; cAtoms+=1){
      //for Cer, Cer1P, SM and HexCer
      for (int cAtoms=20; cAtoms<49; cAtoms+=1){
        int dbs = 0;
        int dbsMax = cDbsCombi.get(cAtoms)+1;
        while (dbs<dbsMax){
        //for sphingoid backbone chain library
//          int totalC = cAtoms;
//          int hAtoms = totalC*2-2*dbs+3;        
//          int oAtoms = 2;

          //for Sphingoid bases
        //int totalC = cAtoms;
        //int hAtoms = totalC*2-2*dbs+3;        
        //for mSphingoid bases
        //int oAtoms = 1;
        //for dSphingoid bases
        //int oAtoms = 2;
        //for tSphingoid bases
        //int oAtoms = r;

          //for Cer
          //int totalC = cAtoms;
          //int hAtoms = totalC*2-2*dbs+1;
          //for mCer
          //int oAtoms = 2;
          //for dCer
          //int oAtoms = 3;
          //for tCer
          //int oAtoms = 4;
          //for qCer
          //int oAtoms = 5;

        //for Cer1P
        //int totalC = cAtoms;
        //int hAtoms = totalC*2-2*dbs+2;
        //for mCer1P
        //int oAtoms = 5;
        //for dCer1P
        //int oAtoms = 6;
        //for tCer1P
        //int oAtoms = 7;

          //for SM
//          int totalC = cAtoms+5;
//          int hAtoms = totalC*2-2*dbs+3;
//          //for mSM
//          //int oAtoms = 5;
//          //for dSM
//          //int oAtoms = 6;
//          //for tSM
//          int oAtoms = 7;

          //for LSM
//          int totalC = cAtoms+5;
//          int hAtoms = totalC*2-2*dbs+5;
//          //for mLSM
//          //int oAtoms = 4;
//          //for dLSM
//          //int oAtoms = 5;
//          //for tLSM
//          int oAtoms = 6;

          //for HexCer
          int totalC = cAtoms+6;
          int hAtoms = totalC*2-2*dbs-1;
//          //for mHexCer
//          //int oAtoms = 7;
//          //for dHexCer
          int oAtoms = 8;
//          //for tHexCer
//          int oAtoms = 9;

          //for LacCer
//          int totalC = cAtoms+12;
//          int hAtoms = totalC*2-2*dbs-3;
//          //for dLacCer
//          int oAtoms = 13;

          //for S1P
//          int totalC = cAtoms;
//          int hAtoms = totalC*2-2*dbs+4;
//          //for mS1P
//          //int oAtoms = 4;
//          //for dS1P
//          //int oAtoms = 5;
//          //for tS1P
//          int oAtoms = 6;
          
          //for sphingoid backbone chain library, sphingoid bases, Cer, HexCer and LacCer
          int pAtoms = 0;
          //for Cer1P, SM, LSM and S1P
          //int pAtoms = 1;
          //for sphingoid backbone chain library, sphingoid bases, Cer, Cer1P, HexCer, LacCer and S1P
          int nAtoms = 1;
          //for SM and LSM
          //int nAtoms = 2;
          int dAtoms = 0;
          double massNeutral = parser.getElementDetails("C").getMonoMass()*((double)totalC)+parser.getElementDetails("H").getMonoMass()*((double)hAtoms)+
              parser.getElementDetails("O").getMonoMass()*((double)oAtoms)+parser.getElementDetails("P").getMonoMass()*((double)pAtoms)+
              parser.getElementDetails("N").getMonoMass()*((double)nAtoms)+parser.getElementDetails("D").getMonoMass()*((double)dAtoms);
          double massCharged = massNeutral+parser.getElementDetails("H").getMonoMass()+parser.getElementDetails("C").getMonoMass()+2d*parser.getElementDetails("O").getMonoMass();
          
          //for positive ion mode
          //protonated
          double massH = massNeutral+parser.getElementDetails("h").getMonoMass();
          //water loss
          double massOH = massNeutral+parser.getElementDetails("h").getMonoMass()-2*parser.getElementDetails("H").getMonoMass()-parser.getElementDetails("O").getMonoMass();
          //sodiated
          double massNa = massNeutral+parser.getElementDetails("na").getMonoMass();

          //for negative ion mode
          //deprotonated
//          double massH = massNeutral-parser.getElementDetails("h").getMonoMass();
//          //water loss
//          double massHCOO = massNeutral-parser.getElementDetails("h").getMonoMass()+2*parser.getElementDetails("H").getMonoMass()+parser.getElementDetails("C").getMonoMass()+2*parser.getElementDetails("O").getMonoMass();
//          //sodiated
//          double massCl = massNeutral+parser.getElementDetails("cl").getMonoMass();

          outRow = resultSheet.createRow(rowCount);
          label = outRow.createCell(0,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(cAtoms);
          label = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
          label.setCellValue(":");       
          label = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(dbs);
          label = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(totalC);
          label = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(hAtoms);
          label = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(oAtoms);
          label = outRow.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(pAtoms);
          label = outRow.createCell(7,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(nAtoms);
          label = outRow.createCell(8,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(dAtoms);
          label = outRow.createCell(9,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(massNeutral);
          
          //positive
          label = outRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(massH);
          label = outRow.createCell(11,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(massOH);
          label = outRow.createCell(12,HSSFCell.CELL_TYPE_NUMERIC);
          label.setCellValue(massNa);
          //negative
//        label = outRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
//        label.setCellValue(massH);
//        label = outRow.createCell(11,HSSFCell.CELL_TYPE_NUMERIC);
//        label.setCellValue(massHCOO);
//        label = outRow.createCell(12,HSSFCell.CELL_TYPE_NUMERIC);
//        label.setCellValue(massCl);

          dbs++;
          rowCount++;
        }
      }
      
      resultWorkbook.write(out);
      out.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }

  }
  
  private void parseLCBList() {
    File lcbFile = new File("C:\\Development\\LipidDataAnalyzer\\fattyAcids\\dLCB.xlsx");
    try {
      LCBLibParser parser = new LCBLibParser(lcbFile);
      parser.parseFile();
    }
    catch (IOException | RulesException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
//    parser.
  }
  
  private void replicateCudaBug() {
    CgChromatogram cgc = new CgChromatogram(800);
    cgc.Value[0][0] = 439.92114f;
    cgc.Value[0][1] = 416666.66f;
    cgc.Value[1][0] = 439.92615f;
    cgc.Value[1][1] = 0.0f;
    cgc.Value[2][0] = 439.93115f;
    cgc.Value[2][1] = 0.0f;
    cgc.Value[3][0] = 439.93616f;
    cgc.Value[3][1] = 0.0f;
    cgc.Value[4][0] = 439.94113f;
    cgc.Value[4][1] = 0.0f;
    cgc.Value[5][0] = 439.94614f;
    cgc.Value[5][1] = 0.0f;
    cgc.Value[6][0] = 439.95114f;
    cgc.Value[6][1] = 0.0f;
    cgc.Value[7][0] = 439.95615f;
    cgc.Value[7][1] = 0.0f;
    cgc.Value[8][0] = 439.96115f;
    cgc.Value[8][1] = 0.0f;
    cgc.Value[9][0] = 439.96616f;
    cgc.Value[9][1] = 0.0f;
    cgc.Value[10][0] = 439.97113f;
    cgc.Value[10][1] = 0.0f;
    cgc.Value[11][0] = 439.97614f;
    cgc.Value[11][1] = 0.0f;
    cgc.Value[12][0] = 439.98114f;
    cgc.Value[12][1] = 0.0f;
    cgc.Value[13][0] = 439.98615f;
    cgc.Value[13][1] = 0.0f;
    cgc.Value[14][0] = 439.99115f;
    cgc.Value[14][1] = 0.0f;
    cgc.Value[15][0] = 439.99615f;
    cgc.Value[15][1] = 0.0f;
    cgc.Value[16][0] = 440.00113f;
    cgc.Value[16][1] = 0.0f;
    cgc.Value[17][0] = 440.00613f;
    cgc.Value[17][1] = 0.0f;
    cgc.Value[18][0] = 440.01114f;
    cgc.Value[18][1] = 0.0f;
    cgc.Value[19][0] = 440.01614f;
    cgc.Value[19][1] = 0.0f;
    cgc.Value[20][0] = 440.02115f;
    cgc.Value[20][1] = 0.0f;
    cgc.Value[21][0] = 440.02615f;
    cgc.Value[21][1] = 0.0f;
    cgc.Value[22][0] = 440.03113f;
    cgc.Value[22][1] = 0.0f;
    cgc.Value[23][0] = 440.03613f;
    cgc.Value[23][1] = 0.0f;
    cgc.Value[24][0] = 440.04114f;
    cgc.Value[24][1] = 83333.336f;
    cgc.Value[25][0] = 440.04614f;
    cgc.Value[25][1] = 0.0f;
    cgc.Value[26][0] = 440.05115f;
    cgc.Value[26][1] = 0.0f;
    cgc.Value[27][0] = 440.05615f;
    cgc.Value[27][1] = 0.0f;
    cgc.Value[28][0] = 440.06116f;
    cgc.Value[28][1] = 0.0f;
    cgc.Value[29][0] = 440.06613f;
    cgc.Value[29][1] = 0.0f;
    cgc.Value[30][0] = 440.07114f;
    cgc.Value[30][1] = 0.0f;
    cgc.Value[31][0] = 440.07614f;
    cgc.Value[31][1] = 0.0f;
    cgc.Value[32][0] = 440.08115f;
    cgc.Value[32][1] = 0.0f;
    cgc.Value[33][0] = 440.08615f;
    cgc.Value[33][1] = 0.0f;
    cgc.Value[34][0] = 440.09116f;
    cgc.Value[34][1] = 0.0f;
    cgc.Value[35][0] = 440.09613f;
    cgc.Value[35][1] = 0.0f;
    cgc.Value[36][0] = 440.10114f;
    cgc.Value[36][1] = 0.0f;
    cgc.Value[37][0] = 440.10614f;
    cgc.Value[37][1] = 0.0f;
    cgc.Value[38][0] = 440.11115f;
    cgc.Value[38][1] = 0.0f;
    cgc.Value[39][0] = 440.11615f;
    cgc.Value[39][1] = 0.0f;
    cgc.Value[40][0] = 440.12115f;
    cgc.Value[40][1] = 0.0f;
    cgc.Value[41][0] = 440.12613f;
    cgc.Value[41][1] = 0.0f;
    cgc.Value[42][0] = 440.13113f;
    cgc.Value[42][1] = 0.0f;
    cgc.Value[43][0] = 440.13614f;
    cgc.Value[43][1] = 0.0f;
    cgc.Value[44][0] = 440.14114f;
    cgc.Value[44][1] = 0.0f;
    cgc.Value[45][0] = 440.14615f;
    cgc.Value[45][1] = 0.0f;
    cgc.Value[46][0] = 440.15115f;
    cgc.Value[46][1] = 0.0f;
    cgc.Value[47][0] = 440.15613f;
    cgc.Value[47][1] = 0.0f;
    cgc.Value[48][0] = 440.16113f;
    cgc.Value[48][1] = 0.0f;
    cgc.Value[49][0] = 440.16614f;
    cgc.Value[49][1] = 0.0f;
    cgc.Value[50][0] = 440.17114f;
    cgc.Value[50][1] = 0.0f;
    cgc.Value[51][0] = 440.17615f;
    cgc.Value[51][1] = 0.0f;
    cgc.Value[52][0] = 440.18115f;
    cgc.Value[52][1] = 0.0f;
    cgc.Value[53][0] = 440.18616f;
    cgc.Value[53][1] = 0.0f;
    cgc.Value[54][0] = 440.19113f;
    cgc.Value[54][1] = 0.0f;
    cgc.Value[55][0] = 440.19614f;
    cgc.Value[55][1] = 0.0f;
    cgc.Value[56][0] = 440.20114f;
    cgc.Value[56][1] = 0.0f;
    cgc.Value[57][0] = 440.20615f;
    cgc.Value[57][1] = 0.0f;
    cgc.Value[58][0] = 440.21115f;
    cgc.Value[58][1] = 0.0f;
    cgc.Value[59][0] = 440.21616f;
    cgc.Value[59][1] = 0.0f;
    cgc.Value[60][0] = 440.22113f;
    cgc.Value[60][1] = 0.0f;
    cgc.Value[61][0] = 440.22614f;
    cgc.Value[61][1] = 0.0f;
    cgc.Value[62][0] = 440.23114f;
    cgc.Value[62][1] = 0.0f;
    cgc.Value[63][0] = 440.23615f;
    cgc.Value[63][1] = 0.0f;
    cgc.Value[64][0] = 440.24115f;
    cgc.Value[64][1] = 0.0f;
    cgc.Value[65][0] = 440.24615f;
    cgc.Value[65][1] = 0.0f;
    cgc.Value[66][0] = 440.25113f;
    cgc.Value[66][1] = 0.0f;
    cgc.Value[67][0] = 440.25613f;
    cgc.Value[67][1] = 0.0f;
    cgc.Value[68][0] = 440.26114f;
    cgc.Value[68][1] = 0.0f;
    cgc.Value[69][0] = 440.26614f;
    cgc.Value[69][1] = 0.0f;
    cgc.Value[70][0] = 440.27115f;
    cgc.Value[70][1] = 0.0f;
    cgc.Value[71][0] = 440.27615f;
    cgc.Value[71][1] = 0.0f;
    cgc.Value[72][0] = 440.28113f;
    cgc.Value[72][1] = 0.0f;
    cgc.Value[73][0] = 440.28613f;
    cgc.Value[73][1] = 0.0f;
    cgc.Value[74][0] = 440.29114f;
    cgc.Value[74][1] = 0.0f;
    cgc.Value[75][0] = 440.29614f;
    cgc.Value[75][1] = 0.0f;
    cgc.Value[76][0] = 440.30115f;
    cgc.Value[76][1] = 0.0f;
    cgc.Value[77][0] = 440.30615f;
    cgc.Value[77][1] = 0.0f;
    cgc.Value[78][0] = 440.31116f;
    cgc.Value[78][1] = 0.0f;
    cgc.Value[79][0] = 440.31613f;
    cgc.Value[79][1] = 0.0f;
    cgc.Value[80][0] = 440.32114f;
    cgc.Value[80][1] = 0.0f;
    cgc.Value[81][0] = 440.32614f;
    cgc.Value[81][1] = 0.0f;
    cgc.Value[82][0] = 440.33115f;
    cgc.Value[82][1] = 0.0f;
    cgc.Value[83][0] = 440.33615f;
    cgc.Value[83][1] = 0.0f;
    cgc.Value[84][0] = 440.34116f;
    cgc.Value[84][1] = 0.0f;
    cgc.Value[85][0] = 440.34613f;
    cgc.Value[85][1] = 0.0f;
    cgc.Value[86][0] = 440.35114f;
    cgc.Value[86][1] = 0.0f;
    cgc.Value[87][0] = 440.35614f;
    cgc.Value[87][1] = 0.0f;
    cgc.Value[88][0] = 440.36115f;
    cgc.Value[88][1] = 0.0f;
    cgc.Value[89][0] = 440.36615f;
    cgc.Value[89][1] = 0.0f;
    cgc.Value[90][0] = 440.37115f;
    cgc.Value[90][1] = 0.0f;
    cgc.Value[91][0] = 440.37613f;
    cgc.Value[91][1] = 0.0f;
    cgc.Value[92][0] = 440.38113f;
    cgc.Value[92][1] = 0.0f;
    cgc.Value[93][0] = 440.38614f;
    cgc.Value[93][1] = 0.0f;
    cgc.Value[94][0] = 440.39114f;
    cgc.Value[94][1] = 0.0f;
    cgc.Value[95][0] = 440.39615f;
    cgc.Value[95][1] = 0.0f;
    cgc.Value[96][0] = 440.40115f;
    cgc.Value[96][1] = 0.0f;
    cgc.Value[97][0] = 440.40613f;
    cgc.Value[97][1] = 0.0f;
    cgc.Value[98][0] = 440.41113f;
    cgc.Value[98][1] = 0.0f;
    cgc.Value[99][0] = 440.41614f;
    cgc.Value[99][1] = 0.0f;
    cgc.Value[100][0] = 440.42114f;
    cgc.Value[100][1] = 0.0f;
    cgc.Value[101][0] = 440.42615f;
    cgc.Value[101][1] = 0.0f;
    cgc.Value[102][0] = 440.43115f;
    cgc.Value[102][1] = 0.0f;
    cgc.Value[103][0] = 440.43616f;
    cgc.Value[103][1] = 0.0f;
    cgc.Value[104][0] = 440.44113f;
    cgc.Value[104][1] = 0.0f;
    cgc.Value[105][0] = 440.44614f;
    cgc.Value[105][1] = 0.0f;
    cgc.Value[106][0] = 440.45114f;
    cgc.Value[106][1] = 0.0f;
    cgc.Value[107][0] = 440.45615f;
    cgc.Value[107][1] = 0.0f;
    cgc.Value[108][0] = 440.46115f;
    cgc.Value[108][1] = 0.0f;
    cgc.Value[109][0] = 440.46616f;
    cgc.Value[109][1] = 0.0f;
    cgc.Value[110][0] = 440.47113f;
    cgc.Value[110][1] = 0.0f;
    cgc.Value[111][0] = 440.47614f;
    cgc.Value[111][1] = 0.0f;
    cgc.Value[112][0] = 440.48114f;
    cgc.Value[112][1] = 0.0f;
    cgc.Value[113][0] = 440.48615f;
    cgc.Value[113][1] = 0.0f;
    cgc.Value[114][0] = 440.49115f;
    cgc.Value[114][1] = 0.0f;
    cgc.Value[115][0] = 440.49615f;
    cgc.Value[115][1] = 0.0f;
    cgc.Value[116][0] = 440.50113f;
    cgc.Value[116][1] = 0.0f;
    cgc.Value[117][0] = 440.50613f;
    cgc.Value[117][1] = 0.0f;
    cgc.Value[118][0] = 440.51114f;
    cgc.Value[118][1] = 0.0f;
    cgc.Value[119][0] = 440.51614f;
    cgc.Value[119][1] = 0.0f;
    cgc.Value[120][0] = 440.52115f;
    cgc.Value[120][1] = 0.0f;
    cgc.Value[121][0] = 440.52615f;
    cgc.Value[121][1] = 0.0f;
    cgc.Value[122][0] = 440.53113f;
    cgc.Value[122][1] = 0.0f;
    cgc.Value[123][0] = 440.53613f;
    cgc.Value[123][1] = 0.0f;
    cgc.Value[124][0] = 440.54114f;
    cgc.Value[124][1] = 0.0f;
    cgc.Value[125][0] = 440.54614f;
    cgc.Value[125][1] = 0.0f;
    cgc.Value[126][0] = 440.55115f;
    cgc.Value[126][1] = 0.0f;
    cgc.Value[127][0] = 440.55615f;
    cgc.Value[127][1] = 0.0f;
    cgc.Value[128][0] = 440.56116f;
    cgc.Value[128][1] = 0.0f;
    cgc.Value[129][0] = 440.56613f;
    cgc.Value[129][1] = 0.0f;
    cgc.Value[130][0] = 440.57114f;
    cgc.Value[130][1] = 0.0f;
    cgc.Value[131][0] = 440.57614f;
    cgc.Value[131][1] = 0.0f;
    cgc.Value[132][0] = 440.58115f;
    cgc.Value[132][1] = 0.0f;
    cgc.Value[133][0] = 440.58615f;
    cgc.Value[133][1] = 0.0f;
    cgc.Value[134][0] = 440.59116f;
    cgc.Value[134][1] = 0.0f;
    cgc.Value[135][0] = 440.59613f;
    cgc.Value[135][1] = 0.0f;
    cgc.Value[136][0] = 440.60114f;
    cgc.Value[136][1] = 0.0f;
    cgc.Value[137][0] = 440.60614f;
    cgc.Value[137][1] = 0.0f;
    cgc.Value[138][0] = 440.61115f;
    cgc.Value[138][1] = 0.0f;
    cgc.Value[139][0] = 440.61615f;
    cgc.Value[139][1] = 0.0f;
    cgc.Value[140][0] = 440.62115f;
    cgc.Value[140][1] = 0.0f;
    cgc.Value[141][0] = 440.62613f;
    cgc.Value[141][1] = 0.0f;
    cgc.Value[142][0] = 440.63113f;
    cgc.Value[142][1] = 0.0f;
    cgc.Value[143][0] = 440.63614f;
    cgc.Value[143][1] = 0.0f;
    cgc.Value[144][0] = 440.64114f;
    cgc.Value[144][1] = 0.0f;
    cgc.Value[145][0] = 440.64615f;
    cgc.Value[145][1] = 0.0f;
    cgc.Value[146][0] = 440.65115f;
    cgc.Value[146][1] = 0.0f;
    cgc.Value[147][0] = 440.65613f;
    cgc.Value[147][1] = 0.0f;
    cgc.Value[148][0] = 440.66113f;
    cgc.Value[148][1] = 0.0f;
    cgc.Value[149][0] = 440.66614f;
    cgc.Value[149][1] = 0.0f;
    cgc.Value[150][0] = 440.67114f;
    cgc.Value[150][1] = 0.0f;
    cgc.Value[151][0] = 440.67615f;
    cgc.Value[151][1] = 0.0f;
    cgc.Value[152][0] = 440.68115f;
    cgc.Value[152][1] = 0.0f;
    cgc.Value[153][0] = 440.68616f;
    cgc.Value[153][1] = 0.0f;
    cgc.Value[154][0] = 440.69113f;
    cgc.Value[154][1] = 0.0f;
    cgc.Value[155][0] = 440.69614f;
    cgc.Value[155][1] = 0.0f;
    cgc.Value[156][0] = 440.70114f;
    cgc.Value[156][1] = 0.0f;
    cgc.Value[157][0] = 440.70615f;
    cgc.Value[157][1] = 0.0f;
    cgc.Value[158][0] = 440.71115f;
    cgc.Value[158][1] = 0.0f;
    cgc.Value[159][0] = 440.71616f;
    cgc.Value[159][1] = 0.0f;
    cgc.Value[160][0] = 440.72113f;
    cgc.Value[160][1] = 0.0f;
    cgc.Value[161][0] = 440.72614f;
    cgc.Value[161][1] = 0.0f;
    cgc.Value[162][0] = 440.73114f;
    cgc.Value[162][1] = 0.0f;
    cgc.Value[163][0] = 440.73615f;
    cgc.Value[163][1] = 0.0f;
    cgc.Value[164][0] = 440.74115f;
    cgc.Value[164][1] = 0.0f;
    cgc.Value[165][0] = 440.74615f;
    cgc.Value[165][1] = 0.0f;
    cgc.Value[166][0] = 440.75113f;
    cgc.Value[166][1] = 0.0f;
    cgc.Value[167][0] = 440.75613f;
    cgc.Value[167][1] = 0.0f;
    cgc.Value[168][0] = 440.76114f;
    cgc.Value[168][1] = 0.0f;
    cgc.Value[169][0] = 440.76614f;
    cgc.Value[169][1] = 0.0f;
    cgc.Value[170][0] = 440.77115f;
    cgc.Value[170][1] = 0.0f;
    cgc.Value[171][0] = 440.77615f;
    cgc.Value[171][1] = 0.0f;
    cgc.Value[172][0] = 440.78113f;
    cgc.Value[172][1] = 0.0f;
    cgc.Value[173][0] = 440.78613f;
    cgc.Value[173][1] = 0.0f;
    cgc.Value[174][0] = 440.79114f;
    cgc.Value[174][1] = 0.0f;
    cgc.Value[175][0] = 440.79614f;
    cgc.Value[175][1] = 0.0f;
    cgc.Value[176][0] = 440.80115f;
    cgc.Value[176][1] = 0.0f;
    cgc.Value[177][0] = 440.80615f;
    cgc.Value[177][1] = 0.0f;
    cgc.Value[178][0] = 440.81116f;
    cgc.Value[178][1] = 0.0f;
    cgc.Value[179][0] = 440.81613f;
    cgc.Value[179][1] = 0.0f;
    cgc.Value[180][0] = 440.82114f;
    cgc.Value[180][1] = 0.0f;
    cgc.Value[181][0] = 440.82614f;
    cgc.Value[181][1] = 0.0f;
    cgc.Value[182][0] = 440.83115f;
    cgc.Value[182][1] = 0.0f;
    cgc.Value[183][0] = 440.83615f;
    cgc.Value[183][1] = 0.0f;
    cgc.Value[184][0] = 440.84116f;
    cgc.Value[184][1] = 0.0f;
    cgc.Value[185][0] = 440.84613f;
    cgc.Value[185][1] = 0.0f;
    cgc.Value[186][0] = 440.85114f;
    cgc.Value[186][1] = 0.0f;
    cgc.Value[187][0] = 440.85614f;
    cgc.Value[187][1] = 0.0f;
    cgc.Value[188][0] = 440.86115f;
    cgc.Value[188][1] = 0.0f;
    cgc.Value[189][0] = 440.86615f;
    cgc.Value[189][1] = 0.0f;
    cgc.Value[190][0] = 440.87115f;
    cgc.Value[190][1] = 0.0f;
    cgc.Value[191][0] = 440.87613f;
    cgc.Value[191][1] = 0.0f;
    cgc.Value[192][0] = 440.88113f;
    cgc.Value[192][1] = 750000.0f;
    cgc.Value[193][0] = 440.88614f;
    cgc.Value[193][1] = 0.0f;
    cgc.Value[194][0] = 440.89114f;
    cgc.Value[194][1] = 0.0f;
    cgc.Value[195][0] = 440.89615f;
    cgc.Value[195][1] = 0.0f;
    cgc.Value[196][0] = 440.90115f;
    cgc.Value[196][1] = 0.0f;
    cgc.Value[197][0] = 440.90613f;
    cgc.Value[197][1] = 0.0f;
    cgc.Value[198][0] = 440.91113f;
    cgc.Value[198][1] = 0.0f;
    cgc.Value[199][0] = 440.91614f;
    cgc.Value[199][1] = 0.0f;
    cgc.Value[200][0] = 440.92114f;
    cgc.Value[200][1] = 0.0f;
    cgc.Value[201][0] = 440.92615f;
    cgc.Value[201][1] = 0.0f;
    cgc.Value[202][0] = 440.93115f;
    cgc.Value[202][1] = 0.0f;
    cgc.Value[203][0] = 440.93616f;
    cgc.Value[203][1] = 0.0f;
    cgc.Value[204][0] = 440.94113f;
    cgc.Value[204][1] = 0.0f;
    cgc.Value[205][0] = 440.94614f;
    cgc.Value[205][1] = 0.0f;
    cgc.Value[206][0] = 440.95114f;
    cgc.Value[206][1] = 0.0f;
    cgc.Value[207][0] = 440.95615f;
    cgc.Value[207][1] = 0.0f;
    cgc.Value[208][0] = 440.96115f;
    cgc.Value[208][1] = 0.0f;
    cgc.Value[209][0] = 440.96616f;
    cgc.Value[209][1] = 0.0f;
    cgc.Value[210][0] = 440.97113f;
    cgc.Value[210][1] = 0.0f;
    cgc.Value[211][0] = 440.97614f;
    cgc.Value[211][1] = 0.0f;
    cgc.Value[212][0] = 440.98114f;
    cgc.Value[212][1] = 0.0f;
    cgc.Value[213][0] = 440.98615f;
    cgc.Value[213][1] = 0.0f;
    cgc.Value[214][0] = 440.99115f;
    cgc.Value[214][1] = 0.0f;
    cgc.Value[215][0] = 440.99615f;
    cgc.Value[215][1] = 0.0f;
    cgc.Value[216][0] = 441.00113f;
    cgc.Value[216][1] = 0.0f;
    cgc.Value[217][0] = 441.00613f;
    cgc.Value[217][1] = 0.0f;
    cgc.Value[218][0] = 441.01114f;
    cgc.Value[218][1] = 0.0f;
    cgc.Value[219][0] = 441.01614f;
    cgc.Value[219][1] = 0.0f;
    cgc.Value[220][0] = 441.02115f;
    cgc.Value[220][1] = 0.0f;
    cgc.Value[221][0] = 441.02615f;
    cgc.Value[221][1] = 0.0f;
    cgc.Value[222][0] = 441.03113f;
    cgc.Value[222][1] = 0.0f;
    cgc.Value[223][0] = 441.03613f;
    cgc.Value[223][1] = 0.0f;
    cgc.Value[224][0] = 441.04114f;
    cgc.Value[224][1] = 0.0f;
    cgc.Value[225][0] = 441.04614f;
    cgc.Value[225][1] = 0.0f;
    cgc.Value[226][0] = 441.05115f;
    cgc.Value[226][1] = 0.0f;
    cgc.Value[227][0] = 441.05615f;
    cgc.Value[227][1] = 0.0f;
    cgc.Value[228][0] = 441.06116f;
    cgc.Value[228][1] = 0.0f;
    cgc.Value[229][0] = 441.06613f;
    cgc.Value[229][1] = 0.0f;
    cgc.Value[230][0] = 441.07114f;
    cgc.Value[230][1] = 0.0f;
    cgc.Value[231][0] = 441.07614f;
    cgc.Value[231][1] = 0.0f;
    cgc.Value[232][0] = 441.08115f;
    cgc.Value[232][1] = 0.0f;
    cgc.Value[233][0] = 441.08615f;
    cgc.Value[233][1] = 0.0f;
    cgc.Value[234][0] = 441.09116f;
    cgc.Value[234][1] = 0.0f;
    cgc.Value[235][0] = 441.09613f;
    cgc.Value[235][1] = 0.0f;
    cgc.Value[236][0] = 441.10114f;
    cgc.Value[236][1] = 0.0f;
    cgc.Value[237][0] = 441.10614f;
    cgc.Value[237][1] = 0.0f;
    cgc.Value[238][0] = 441.11115f;
    cgc.Value[238][1] = 0.0f;
    cgc.Value[239][0] = 441.11615f;
    cgc.Value[239][1] = 0.0f;
    cgc.Value[240][0] = 441.12115f;
    cgc.Value[240][1] = 166666.67f;
    cgc.Value[241][0] = 441.12613f;
    cgc.Value[241][1] = 0.0f;
    cgc.Value[242][0] = 441.13113f;
    cgc.Value[242][1] = 0.0f;
    cgc.Value[243][0] = 441.13614f;
    cgc.Value[243][1] = 0.0f;
    cgc.Value[244][0] = 441.14114f;
    cgc.Value[244][1] = 0.0f;
    cgc.Value[245][0] = 441.14615f;
    cgc.Value[245][1] = 0.0f;
    cgc.Value[246][0] = 441.15115f;
    cgc.Value[246][1] = 0.0f;
    cgc.Value[247][0] = 441.15613f;
    cgc.Value[247][1] = 0.0f;
    cgc.Value[248][0] = 441.16113f;
    cgc.Value[248][1] = 0.0f;
    cgc.Value[249][0] = 441.16614f;
    cgc.Value[249][1] = 0.0f;
    cgc.Value[250][0] = 441.17114f;
    cgc.Value[250][1] = 0.0f;
    cgc.Value[251][0] = 441.17615f;
    cgc.Value[251][1] = 0.0f;
    cgc.Value[252][0] = 441.18115f;
    cgc.Value[252][1] = 0.0f;
    cgc.Value[253][0] = 441.18616f;
    cgc.Value[253][1] = 0.0f;
    cgc.Value[254][0] = 441.19113f;
    cgc.Value[254][1] = 0.0f;
    cgc.Value[255][0] = 441.19614f;
    cgc.Value[255][1] = 0.0f;
    cgc.Value[256][0] = 441.20114f;
    cgc.Value[256][1] = 0.0f;
    cgc.Value[257][0] = 441.20615f;
    cgc.Value[257][1] = 0.0f;
    cgc.Value[258][0] = 441.21115f;
    cgc.Value[258][1] = 0.0f;
    cgc.Value[259][0] = 441.21616f;
    cgc.Value[259][1] = 0.0f;
    cgc.Value[260][0] = 441.22113f;
    cgc.Value[260][1] = 0.0f;
    cgc.Value[261][0] = 441.22614f;
    cgc.Value[261][1] = 0.0f;
    cgc.Value[262][0] = 441.23114f;
    cgc.Value[262][1] = 0.0f;
    cgc.Value[263][0] = 441.23615f;
    cgc.Value[263][1] = 0.0f;
    cgc.Value[264][0] = 441.24115f;
    cgc.Value[264][1] = 0.0f;
    cgc.Value[265][0] = 441.24615f;
    cgc.Value[265][1] = 0.0f;
    cgc.Value[266][0] = 441.25113f;
    cgc.Value[266][1] = 0.0f;
    cgc.Value[267][0] = 441.25613f;
    cgc.Value[267][1] = 0.0f;
    cgc.Value[268][0] = 441.26114f;
    cgc.Value[268][1] = 0.0f;
    cgc.Value[269][0] = 441.26614f;
    cgc.Value[269][1] = 0.0f;
    cgc.Value[270][0] = 441.27115f;
    cgc.Value[270][1] = 0.0f;
    cgc.Value[271][0] = 441.27615f;
    cgc.Value[271][1] = 0.0f;
    cgc.Value[272][0] = 441.28113f;
    cgc.Value[272][1] = 0.0f;
    cgc.Value[273][0] = 441.28613f;
    cgc.Value[273][1] = 0.0f;
    cgc.Value[274][0] = 441.29114f;
    cgc.Value[274][1] = 0.0f;
    cgc.Value[275][0] = 441.29614f;
    cgc.Value[275][1] = 0.0f;
    cgc.Value[276][0] = 441.30115f;
    cgc.Value[276][1] = 0.0f;
    cgc.Value[277][0] = 441.30615f;
    cgc.Value[277][1] = 0.0f;
    cgc.Value[278][0] = 441.31116f;
    cgc.Value[278][1] = 0.0f;
    cgc.Value[279][0] = 441.31613f;
    cgc.Value[279][1] = 0.0f;
    cgc.Value[280][0] = 441.32114f;
    cgc.Value[280][1] = 0.0f;
    cgc.Value[281][0] = 441.32614f;
    cgc.Value[281][1] = 0.0f;
    cgc.Value[282][0] = 441.33115f;
    cgc.Value[282][1] = 0.0f;
    cgc.Value[283][0] = 441.33615f;
    cgc.Value[283][1] = 0.0f;
    cgc.Value[284][0] = 441.34116f;
    cgc.Value[284][1] = 0.0f;
    cgc.Value[285][0] = 441.34613f;
    cgc.Value[285][1] = 0.0f;
    cgc.Value[286][0] = 441.35114f;
    cgc.Value[286][1] = 0.0f;
    cgc.Value[287][0] = 441.35614f;
    cgc.Value[287][1] = 0.0f;
    cgc.Value[288][0] = 441.36115f;
    cgc.Value[288][1] = 0.0f;
    cgc.Value[289][0] = 441.36615f;
    cgc.Value[289][1] = 0.0f;
    cgc.Value[290][0] = 441.37115f;
    cgc.Value[290][1] = 0.0f;
    cgc.Value[291][0] = 441.37613f;
    cgc.Value[291][1] = 0.0f;
    cgc.Value[292][0] = 441.38113f;
    cgc.Value[292][1] = 0.0f;
    cgc.Value[293][0] = 441.38614f;
    cgc.Value[293][1] = 0.0f;
    cgc.Value[294][0] = 441.39114f;
    cgc.Value[294][1] = 0.0f;
    cgc.Value[295][0] = 441.39615f;
    cgc.Value[295][1] = 0.0f;
    cgc.Value[296][0] = 441.40115f;
    cgc.Value[296][1] = 0.0f;
    cgc.Value[297][0] = 441.40613f;
    cgc.Value[297][1] = 0.0f;
    cgc.Value[298][0] = 441.41113f;
    cgc.Value[298][1] = 0.0f;
    cgc.Value[299][0] = 441.41614f;
    cgc.Value[299][1] = 0.0f;
    cgc.Value[300][0] = 441.42114f;
    cgc.Value[300][1] = 0.0f;
    cgc.Value[301][0] = 441.42615f;
    cgc.Value[301][1] = 0.0f;
    cgc.Value[302][0] = 441.43115f;
    cgc.Value[302][1] = 0.0f;
    cgc.Value[303][0] = 441.43616f;
    cgc.Value[303][1] = 0.0f;
    cgc.Value[304][0] = 441.44113f;
    cgc.Value[304][1] = 0.0f;
    cgc.Value[305][0] = 441.44614f;
    cgc.Value[305][1] = 0.0f;
    cgc.Value[306][0] = 441.45114f;
    cgc.Value[306][1] = 0.0f;
    cgc.Value[307][0] = 441.45615f;
    cgc.Value[307][1] = 0.0f;
    cgc.Value[308][0] = 441.46115f;
    cgc.Value[308][1] = 0.0f;
    cgc.Value[309][0] = 441.46616f;
    cgc.Value[309][1] = 0.0f;
    cgc.Value[310][0] = 441.47113f;
    cgc.Value[310][1] = 0.0f;
    cgc.Value[311][0] = 441.47614f;
    cgc.Value[311][1] = 0.0f;
    cgc.Value[312][0] = 441.48114f;
    cgc.Value[312][1] = 0.0f;
    cgc.Value[313][0] = 441.48615f;
    cgc.Value[313][1] = 0.0f;
    cgc.Value[314][0] = 441.49115f;
    cgc.Value[314][1] = 0.0f;
    cgc.Value[315][0] = 441.49615f;
    cgc.Value[315][1] = 0.0f;
    cgc.Value[316][0] = 441.50113f;
    cgc.Value[316][1] = 0.0f;
    cgc.Value[317][0] = 441.50613f;
    cgc.Value[317][1] = 0.0f;
    cgc.Value[318][0] = 441.51114f;
    cgc.Value[318][1] = 0.0f;
    cgc.Value[319][0] = 441.51614f;
    cgc.Value[319][1] = 0.0f;
    cgc.Value[320][0] = 441.52115f;
    cgc.Value[320][1] = 0.0f;
    cgc.Value[321][0] = 441.52615f;
    cgc.Value[321][1] = 0.0f;
    cgc.Value[322][0] = 441.53113f;
    cgc.Value[322][1] = 0.0f;
    cgc.Value[323][0] = 441.53613f;
    cgc.Value[323][1] = 0.0f;
    cgc.Value[324][0] = 441.54114f;
    cgc.Value[324][1] = 0.0f;
    cgc.Value[325][0] = 441.54614f;
    cgc.Value[325][1] = 0.0f;
    cgc.Value[326][0] = 441.55115f;
    cgc.Value[326][1] = 0.0f;
    cgc.Value[327][0] = 441.55615f;
    cgc.Value[327][1] = 0.0f;
    cgc.Value[328][0] = 441.56116f;
    cgc.Value[328][1] = 0.0f;
    cgc.Value[329][0] = 441.56613f;
    cgc.Value[329][1] = 0.0f;
    cgc.Value[330][0] = 441.57114f;
    cgc.Value[330][1] = 0.0f;
    cgc.Value[331][0] = 441.57614f;
    cgc.Value[331][1] = 0.0f;
    cgc.Value[332][0] = 441.58115f;
    cgc.Value[332][1] = 0.0f;
    cgc.Value[333][0] = 441.58615f;
    cgc.Value[333][1] = 0.0f;
    cgc.Value[334][0] = 441.59116f;
    cgc.Value[334][1] = 0.0f;
    cgc.Value[335][0] = 441.59613f;
    cgc.Value[335][1] = 0.0f;
    cgc.Value[336][0] = 441.60114f;
    cgc.Value[336][1] = 0.0f;
    cgc.Value[337][0] = 441.60614f;
    cgc.Value[337][1] = 0.0f;
    cgc.Value[338][0] = 441.61115f;
    cgc.Value[338][1] = 0.0f;
    cgc.Value[339][0] = 441.61615f;
    cgc.Value[339][1] = 0.0f;
    cgc.Value[340][0] = 441.62115f;
    cgc.Value[340][1] = 0.0f;
    cgc.Value[341][0] = 441.62613f;
    cgc.Value[341][1] = 0.0f;
    cgc.Value[342][0] = 441.63113f;
    cgc.Value[342][1] = 0.0f;
    cgc.Value[343][0] = 441.63614f;
    cgc.Value[343][1] = 0.0f;
    cgc.Value[344][0] = 441.64114f;
    cgc.Value[344][1] = 0.0f;
    cgc.Value[345][0] = 441.64615f;
    cgc.Value[345][1] = 0.0f;
    cgc.Value[346][0] = 441.65115f;
    cgc.Value[346][1] = 0.0f;
    cgc.Value[347][0] = 441.65613f;
    cgc.Value[347][1] = 0.0f;
    cgc.Value[348][0] = 441.66113f;
    cgc.Value[348][1] = 0.0f;
    cgc.Value[349][0] = 441.66614f;
    cgc.Value[349][1] = 0.0f;
    cgc.Value[350][0] = 441.67114f;
    cgc.Value[350][1] = 0.0f;
    cgc.Value[351][0] = 441.67615f;
    cgc.Value[351][1] = 0.0f;
    cgc.Value[352][0] = 441.68115f;
    cgc.Value[352][1] = 0.0f;
    cgc.Value[353][0] = 441.68616f;
    cgc.Value[353][1] = 0.0f;
    cgc.Value[354][0] = 441.69113f;
    cgc.Value[354][1] = 0.0f;
    cgc.Value[355][0] = 441.69614f;
    cgc.Value[355][1] = 0.0f;
    cgc.Value[356][0] = 441.70114f;
    cgc.Value[356][1] = 0.0f;
    cgc.Value[357][0] = 441.70615f;
    cgc.Value[357][1] = 0.0f;
    cgc.Value[358][0] = 441.71115f;
    cgc.Value[358][1] = 0.0f;
    cgc.Value[359][0] = 441.71616f;
    cgc.Value[359][1] = 0.0f;
    cgc.Value[360][0] = 441.72113f;
    cgc.Value[360][1] = 583333.3f;
    cgc.Value[361][0] = 441.72614f;
    cgc.Value[361][1] = 0.0f;
    cgc.Value[362][0] = 441.73114f;
    cgc.Value[362][1] = 0.0f;
    cgc.Value[363][0] = 441.73615f;
    cgc.Value[363][1] = 0.0f;
    cgc.Value[364][0] = 441.74115f;
    cgc.Value[364][1] = 0.0f;
    cgc.Value[365][0] = 441.74615f;
    cgc.Value[365][1] = 0.0f;
    cgc.Value[366][0] = 441.75113f;
    cgc.Value[366][1] = 0.0f;
    cgc.Value[367][0] = 441.75613f;
    cgc.Value[367][1] = 0.0f;
    cgc.Value[368][0] = 441.76114f;
    cgc.Value[368][1] = 0.0f;
    cgc.Value[369][0] = 441.76614f;
    cgc.Value[369][1] = 0.0f;
    cgc.Value[370][0] = 441.77115f;
    cgc.Value[370][1] = 0.0f;
    cgc.Value[371][0] = 441.77615f;
    cgc.Value[371][1] = 0.0f;
    cgc.Value[372][0] = 441.78113f;
    cgc.Value[372][1] = 0.0f;
    cgc.Value[373][0] = 441.78613f;
    cgc.Value[373][1] = 0.0f;
    cgc.Value[374][0] = 441.79114f;
    cgc.Value[374][1] = 0.0f;
    cgc.Value[375][0] = 441.79614f;
    cgc.Value[375][1] = 0.0f;
    cgc.Value[376][0] = 441.80115f;
    cgc.Value[376][1] = 0.0f;
    cgc.Value[377][0] = 441.80615f;
    cgc.Value[377][1] = 0.0f;
    cgc.Value[378][0] = 441.81116f;
    cgc.Value[378][1] = 0.0f;
    cgc.Value[379][0] = 441.81613f;
    cgc.Value[379][1] = 0.0f;
    cgc.Value[380][0] = 441.82114f;
    cgc.Value[380][1] = 0.0f;
    cgc.Value[381][0] = 441.82614f;
    cgc.Value[381][1] = 0.0f;
    cgc.Value[382][0] = 441.83115f;
    cgc.Value[382][1] = 0.0f;
    cgc.Value[383][0] = 441.83615f;
    cgc.Value[383][1] = 0.0f;
    cgc.Value[384][0] = 441.84116f;
    cgc.Value[384][1] = 250000.0f;
    cgc.Value[385][0] = 441.84613f;
    cgc.Value[385][1] = 0.0f;
    cgc.Value[386][0] = 441.85114f;
    cgc.Value[386][1] = 0.0f;
    cgc.Value[387][0] = 441.85614f;
    cgc.Value[387][1] = 0.0f;
    cgc.Value[388][0] = 441.86115f;
    cgc.Value[388][1] = 0.0f;
    cgc.Value[389][0] = 441.86615f;
    cgc.Value[389][1] = 0.0f;
    cgc.Value[390][0] = 441.87115f;
    cgc.Value[390][1] = 0.0f;
    cgc.Value[391][0] = 441.87613f;
    cgc.Value[391][1] = 0.0f;
    cgc.Value[392][0] = 441.88113f;
    cgc.Value[392][1] = 0.0f;
    cgc.Value[393][0] = 441.88614f;
    cgc.Value[393][1] = 0.0f;
    cgc.Value[394][0] = 441.89114f;
    cgc.Value[394][1] = 0.0f;
    cgc.Value[395][0] = 441.89615f;
    cgc.Value[395][1] = 0.0f;
    cgc.Value[396][0] = 441.90115f;
    cgc.Value[396][1] = 0.0f;
    cgc.Value[397][0] = 441.90613f;
    cgc.Value[397][1] = 0.0f;
    cgc.Value[398][0] = 441.91113f;
    cgc.Value[398][1] = 0.0f;
    cgc.Value[399][0] = 441.91614f;
    cgc.Value[399][1] = 0.0f;
    cgc.Value[400][0] = 441.92114f;
    cgc.Value[400][1] = 0.0f;
    cgc.Value[401][0] = 441.92615f;
    cgc.Value[401][1] = 0.0f;
    cgc.Value[402][0] = 441.93115f;
    cgc.Value[402][1] = 0.0f;
    cgc.Value[403][0] = 441.93616f;
    cgc.Value[403][1] = 0.0f;
    cgc.Value[404][0] = 441.94113f;
    cgc.Value[404][1] = 0.0f;
    cgc.Value[405][0] = 441.94614f;
    cgc.Value[405][1] = 0.0f;
    cgc.Value[406][0] = 441.95114f;
    cgc.Value[406][1] = 0.0f;
    cgc.Value[407][0] = 441.95615f;
    cgc.Value[407][1] = 0.0f;
    cgc.Value[408][0] = 441.96115f;
    cgc.Value[408][1] = 0.0f;
    cgc.Value[409][0] = 441.96616f;
    cgc.Value[409][1] = 0.0f;
    cgc.Value[410][0] = 441.97113f;
    cgc.Value[410][1] = 0.0f;
    cgc.Value[411][0] = 441.97614f;
    cgc.Value[411][1] = 0.0f;
    cgc.Value[412][0] = 441.98114f;
    cgc.Value[412][1] = 0.0f;
    cgc.Value[413][0] = 441.98615f;
    cgc.Value[413][1] = 0.0f;
    cgc.Value[414][0] = 441.99115f;
    cgc.Value[414][1] = 0.0f;
    cgc.Value[415][0] = 441.99615f;
    cgc.Value[415][1] = 0.0f;
    cgc.Value[416][0] = 442.00113f;
    cgc.Value[416][1] = 0.0f;
    cgc.Value[417][0] = 442.00613f;
    cgc.Value[417][1] = 0.0f;
    cgc.Value[418][0] = 442.01114f;
    cgc.Value[418][1] = 0.0f;
    cgc.Value[419][0] = 442.01614f;
    cgc.Value[419][1] = 0.0f;
    cgc.Value[420][0] = 442.02115f;
    cgc.Value[420][1] = 0.0f;
    cgc.Value[421][0] = 442.02615f;
    cgc.Value[421][1] = 0.0f;
    cgc.Value[422][0] = 442.03113f;
    cgc.Value[422][1] = 0.0f;
    cgc.Value[423][0] = 442.03613f;
    cgc.Value[423][1] = 0.0f;
    cgc.Value[424][0] = 442.04114f;
    cgc.Value[424][1] = 0.0f;
    cgc.Value[425][0] = 442.04614f;
    cgc.Value[425][1] = 0.0f;
    cgc.Value[426][0] = 442.05115f;
    cgc.Value[426][1] = 0.0f;
    cgc.Value[427][0] = 442.05615f;
    cgc.Value[427][1] = 0.0f;
    cgc.Value[428][0] = 442.06116f;
    cgc.Value[428][1] = 0.0f;
    cgc.Value[429][0] = 442.06613f;
    cgc.Value[429][1] = 0.0f;
    cgc.Value[430][0] = 442.07114f;
    cgc.Value[430][1] = 0.0f;
    cgc.Value[431][0] = 442.07614f;
    cgc.Value[431][1] = 83333.336f;
    cgc.Value[432][0] = 442.08115f;
    cgc.Value[432][1] = 0.0f;
    cgc.Value[433][0] = 442.08615f;
    cgc.Value[433][1] = 0.0f;
    cgc.Value[434][0] = 442.09116f;
    cgc.Value[434][1] = 0.0f;
    cgc.Value[435][0] = 442.09613f;
    cgc.Value[435][1] = 0.0f;
    cgc.Value[436][0] = 442.10114f;
    cgc.Value[436][1] = 0.0f;
    cgc.Value[437][0] = 442.10614f;
    cgc.Value[437][1] = 0.0f;
    cgc.Value[438][0] = 442.11115f;
    cgc.Value[438][1] = 0.0f;
    cgc.Value[439][0] = 442.11615f;
    cgc.Value[439][1] = 0.0f;
    cgc.Value[440][0] = 442.12115f;
    cgc.Value[440][1] = 0.0f;
    cgc.Value[441][0] = 442.12613f;
    cgc.Value[441][1] = 0.0f;
    cgc.Value[442][0] = 442.13113f;
    cgc.Value[442][1] = 0.0f;
    cgc.Value[443][0] = 442.13614f;
    cgc.Value[443][1] = 0.0f;
    cgc.Value[444][0] = 442.14114f;
    cgc.Value[444][1] = 0.0f;
    cgc.Value[445][0] = 442.14615f;
    cgc.Value[445][1] = 0.0f;
    cgc.Value[446][0] = 442.15115f;
    cgc.Value[446][1] = 0.0f;
    cgc.Value[447][0] = 442.15613f;
    cgc.Value[447][1] = 0.0f;
    cgc.Value[448][0] = 442.16113f;
    cgc.Value[448][1] = 0.0f;
    cgc.Value[449][0] = 442.16614f;
    cgc.Value[449][1] = 0.0f;
    cgc.Value[450][0] = 442.17114f;
    cgc.Value[450][1] = 0.0f;
    cgc.Value[451][0] = 442.17615f;
    cgc.Value[451][1] = 0.0f;
    cgc.Value[452][0] = 442.18115f;
    cgc.Value[452][1] = 0.0f;
    cgc.Value[453][0] = 442.18616f;
    cgc.Value[453][1] = 0.0f;
    cgc.Value[454][0] = 442.19113f;
    cgc.Value[454][1] = 0.0f;
    cgc.Value[455][0] = 442.19614f;
    cgc.Value[455][1] = 0.0f;
    cgc.Value[456][0] = 442.20114f;
    cgc.Value[456][1] = 0.0f;
    cgc.Value[457][0] = 442.20615f;
    cgc.Value[457][1] = 0.0f;
    cgc.Value[458][0] = 442.21115f;
    cgc.Value[458][1] = 0.0f;
    cgc.Value[459][0] = 442.21616f;
    cgc.Value[459][1] = 0.0f;
    cgc.Value[460][0] = 442.22113f;
    cgc.Value[460][1] = 0.0f;
    cgc.Value[461][0] = 442.22614f;
    cgc.Value[461][1] = 0.0f;
    cgc.Value[462][0] = 442.23114f;
    cgc.Value[462][1] = 0.0f;
    cgc.Value[463][0] = 442.23615f;
    cgc.Value[463][1] = 0.0f;
    cgc.Value[464][0] = 442.24115f;
    cgc.Value[464][1] = 0.0f;
    cgc.Value[465][0] = 442.24615f;
    cgc.Value[465][1] = 0.0f;
    cgc.Value[466][0] = 442.25113f;
    cgc.Value[466][1] = 0.0f;
    cgc.Value[467][0] = 442.25613f;
    cgc.Value[467][1] = 0.0f;
    cgc.Value[468][0] = 442.26114f;
    cgc.Value[468][1] = 0.0f;
    cgc.Value[469][0] = 442.26614f;
    cgc.Value[469][1] = 0.0f;
    cgc.Value[470][0] = 442.27115f;
    cgc.Value[470][1] = 0.0f;
    cgc.Value[471][0] = 442.27615f;
    cgc.Value[471][1] = 0.0f;
    cgc.Value[472][0] = 442.28113f;
    cgc.Value[472][1] = 0.0f;
    cgc.Value[473][0] = 442.28613f;
    cgc.Value[473][1] = 0.0f;
    cgc.Value[474][0] = 442.29114f;
    cgc.Value[474][1] = 0.0f;
    cgc.Value[475][0] = 442.29614f;
    cgc.Value[475][1] = 0.0f;
    cgc.Value[476][0] = 442.30115f;
    cgc.Value[476][1] = 0.0f;
    cgc.Value[477][0] = 442.30615f;
    cgc.Value[477][1] = 0.0f;
    cgc.Value[478][0] = 442.31116f;
    cgc.Value[478][1] = 0.0f;
    cgc.Value[479][0] = 442.31613f;
    cgc.Value[479][1] = 0.0f;
    cgc.Value[480][0] = 442.32114f;
    cgc.Value[480][1] = 83333.336f;
    cgc.Value[481][0] = 442.32614f;
    cgc.Value[481][1] = 0.0f;
    cgc.Value[482][0] = 442.33115f;
    cgc.Value[482][1] = 0.0f;
    cgc.Value[483][0] = 442.33615f;
    cgc.Value[483][1] = 0.0f;
    cgc.Value[484][0] = 442.34116f;
    cgc.Value[484][1] = 0.0f;
    cgc.Value[485][0] = 442.34613f;
    cgc.Value[485][1] = 0.0f;
    cgc.Value[486][0] = 442.35114f;
    cgc.Value[486][1] = 0.0f;
    cgc.Value[487][0] = 442.35614f;
    cgc.Value[487][1] = 0.0f;
    cgc.Value[488][0] = 442.36115f;
    cgc.Value[488][1] = 0.0f;
    cgc.Value[489][0] = 442.36615f;
    cgc.Value[489][1] = 0.0f;
    cgc.Value[490][0] = 442.37115f;
    cgc.Value[490][1] = 0.0f;
    cgc.Value[491][0] = 442.37613f;
    cgc.Value[491][1] = 0.0f;
    cgc.Value[492][0] = 442.38113f;
    cgc.Value[492][1] = 0.0f;
    cgc.Value[493][0] = 442.38614f;
    cgc.Value[493][1] = 0.0f;
    cgc.Value[494][0] = 442.39114f;
    cgc.Value[494][1] = 0.0f;
    cgc.Value[495][0] = 442.39615f;
    cgc.Value[495][1] = 0.0f;
    cgc.Value[496][0] = 442.40115f;
    cgc.Value[496][1] = 0.0f;
    cgc.Value[497][0] = 442.40613f;
    cgc.Value[497][1] = 0.0f;
    cgc.Value[498][0] = 442.41113f;
    cgc.Value[498][1] = 0.0f;
    cgc.Value[499][0] = 442.41614f;
    cgc.Value[499][1] = 0.0f;
    cgc.Value[500][0] = 442.42114f;
    cgc.Value[500][1] = 0.0f;
    cgc.Value[501][0] = 442.42615f;
    cgc.Value[501][1] = 0.0f;
    cgc.Value[502][0] = 442.43115f;
    cgc.Value[502][1] = 0.0f;
    cgc.Value[503][0] = 442.43616f;
    cgc.Value[503][1] = 0.0f;
    cgc.Value[504][0] = 442.44113f;
    cgc.Value[504][1] = 0.0f;
    cgc.Value[505][0] = 442.44614f;
    cgc.Value[505][1] = 0.0f;
    cgc.Value[506][0] = 442.45114f;
    cgc.Value[506][1] = 0.0f;
    cgc.Value[507][0] = 442.45615f;
    cgc.Value[507][1] = 0.0f;
    cgc.Value[508][0] = 442.46115f;
    cgc.Value[508][1] = 0.0f;
    cgc.Value[509][0] = 442.46616f;
    cgc.Value[509][1] = 0.0f;
    cgc.Value[510][0] = 442.47113f;
    cgc.Value[510][1] = 0.0f;
    cgc.Value[511][0] = 442.47614f;
    cgc.Value[511][1] = 0.0f;
    cgc.Value[512][0] = 442.48114f;
    cgc.Value[512][1] = 0.0f;
    cgc.Value[513][0] = 442.48615f;
    cgc.Value[513][1] = 0.0f;
    cgc.Value[514][0] = 442.49115f;
    cgc.Value[514][1] = 0.0f;
    cgc.Value[515][0] = 442.49615f;
    cgc.Value[515][1] = 0.0f;
    cgc.Value[516][0] = 442.50113f;
    cgc.Value[516][1] = 0.0f;
    cgc.Value[517][0] = 442.50613f;
    cgc.Value[517][1] = 0.0f;
    cgc.Value[518][0] = 442.51114f;
    cgc.Value[518][1] = 0.0f;
    cgc.Value[519][0] = 442.51614f;
    cgc.Value[519][1] = 0.0f;
    cgc.Value[520][0] = 442.52115f;
    cgc.Value[520][1] = 0.0f;
    cgc.Value[521][0] = 442.52615f;
    cgc.Value[521][1] = 0.0f;
    cgc.Value[522][0] = 442.53113f;
    cgc.Value[522][1] = 0.0f;
    cgc.Value[523][0] = 442.53613f;
    cgc.Value[523][1] = 0.0f;
    cgc.Value[524][0] = 442.54114f;
    cgc.Value[524][1] = 0.0f;
    cgc.Value[525][0] = 442.54614f;
    cgc.Value[525][1] = 0.0f;
    cgc.Value[526][0] = 442.55115f;
    cgc.Value[526][1] = 0.0f;
    cgc.Value[527][0] = 442.55615f;
    cgc.Value[527][1] = 0.0f;
    cgc.Value[528][0] = 442.56116f;
    cgc.Value[528][1] = 166666.67f;
    cgc.Value[529][0] = 442.56613f;
    cgc.Value[529][1] = 0.0f;
    cgc.Value[530][0] = 442.57114f;
    cgc.Value[530][1] = 0.0f;
    cgc.Value[531][0] = 442.57614f;
    cgc.Value[531][1] = 0.0f;
    cgc.Value[532][0] = 442.58115f;
    cgc.Value[532][1] = 0.0f;
    cgc.Value[533][0] = 442.58615f;
    cgc.Value[533][1] = 0.0f;
    cgc.Value[534][0] = 442.59116f;
    cgc.Value[534][1] = 0.0f;
    cgc.Value[535][0] = 442.59613f;
    cgc.Value[535][1] = 0.0f;
    cgc.Value[536][0] = 442.60114f;
    cgc.Value[536][1] = 0.0f;
    cgc.Value[537][0] = 442.60614f;
    cgc.Value[537][1] = 0.0f;
    cgc.Value[538][0] = 442.61115f;
    cgc.Value[538][1] = 0.0f;
    cgc.Value[539][0] = 442.61615f;
    cgc.Value[539][1] = 0.0f;
    cgc.Value[540][0] = 442.62115f;
    cgc.Value[540][1] = 0.0f;
    cgc.Value[541][0] = 442.62613f;
    cgc.Value[541][1] = 0.0f;
    cgc.Value[542][0] = 442.63113f;
    cgc.Value[542][1] = 0.0f;
    cgc.Value[543][0] = 442.63614f;
    cgc.Value[543][1] = 0.0f;
    cgc.Value[544][0] = 442.64114f;
    cgc.Value[544][1] = 0.0f;
    cgc.Value[545][0] = 442.64615f;
    cgc.Value[545][1] = 0.0f;
    cgc.Value[546][0] = 442.65115f;
    cgc.Value[546][1] = 0.0f;
    cgc.Value[547][0] = 442.65613f;
    cgc.Value[547][1] = 0.0f;
    cgc.Value[548][0] = 442.66113f;
    cgc.Value[548][1] = 0.0f;
    cgc.Value[549][0] = 442.66614f;
    cgc.Value[549][1] = 0.0f;
    cgc.Value[550][0] = 442.67114f;
    cgc.Value[550][1] = 0.0f;
    cgc.Value[551][0] = 442.67615f;
    cgc.Value[551][1] = 0.0f;
    cgc.Value[552][0] = 442.68115f;
    cgc.Value[552][1] = 166666.67f;
    cgc.Value[553][0] = 442.68616f;
    cgc.Value[553][1] = 0.0f;
    cgc.Value[554][0] = 442.69113f;
    cgc.Value[554][1] = 0.0f;
    cgc.Value[555][0] = 442.69614f;
    cgc.Value[555][1] = 0.0f;
    cgc.Value[556][0] = 442.70114f;
    cgc.Value[556][1] = 0.0f;
    cgc.Value[557][0] = 442.70615f;
    cgc.Value[557][1] = 0.0f;
    cgc.Value[558][0] = 442.71115f;
    cgc.Value[558][1] = 0.0f;
    cgc.Value[559][0] = 442.71616f;
    cgc.Value[559][1] = 0.0f;
    cgc.Value[560][0] = 442.72113f;
    cgc.Value[560][1] = 0.0f;
    cgc.Value[561][0] = 442.72614f;
    cgc.Value[561][1] = 0.0f;
    cgc.Value[562][0] = 442.73114f;
    cgc.Value[562][1] = 0.0f;
    cgc.Value[563][0] = 442.73615f;
    cgc.Value[563][1] = 0.0f;
    cgc.Value[564][0] = 442.74115f;
    cgc.Value[564][1] = 0.0f;
    cgc.Value[565][0] = 442.74615f;
    cgc.Value[565][1] = 0.0f;
    cgc.Value[566][0] = 442.75113f;
    cgc.Value[566][1] = 0.0f;
    cgc.Value[567][0] = 442.75613f;
    cgc.Value[567][1] = 0.0f;
    cgc.Value[568][0] = 442.76114f;
    cgc.Value[568][1] = 0.0f;
    cgc.Value[569][0] = 442.76614f;
    cgc.Value[569][1] = 0.0f;
    cgc.Value[570][0] = 442.77115f;
    cgc.Value[570][1] = 0.0f;
    cgc.Value[571][0] = 442.77615f;
    cgc.Value[571][1] = 0.0f;
    cgc.Value[572][0] = 442.78113f;
    cgc.Value[572][1] = 0.0f;
    cgc.Value[573][0] = 442.78613f;
    cgc.Value[573][1] = 0.0f;
    cgc.Value[574][0] = 442.79114f;
    cgc.Value[574][1] = 0.0f;
    cgc.Value[575][0] = 442.79614f;
    cgc.Value[575][1] = 250000.0f;
    cgc.Value[576][0] = 442.80115f;
    cgc.Value[576][1] = 0.0f;
    cgc.Value[577][0] = 442.80615f;
    cgc.Value[577][1] = 0.0f;
    cgc.Value[578][0] = 442.81116f;
    cgc.Value[578][1] = 0.0f;
    cgc.Value[579][0] = 442.81613f;
    cgc.Value[579][1] = 0.0f;
    cgc.Value[580][0] = 442.82114f;
    cgc.Value[580][1] = 0.0f;
    cgc.Value[581][0] = 442.82614f;
    cgc.Value[581][1] = 0.0f;
    cgc.Value[582][0] = 442.83115f;
    cgc.Value[582][1] = 0.0f;
    cgc.Value[583][0] = 442.83615f;
    cgc.Value[583][1] = 0.0f;
    cgc.Value[584][0] = 442.84116f;
    cgc.Value[584][1] = 0.0f;
    cgc.Value[585][0] = 442.84613f;
    cgc.Value[585][1] = 0.0f;
    cgc.Value[586][0] = 442.85114f;
    cgc.Value[586][1] = 0.0f;
    cgc.Value[587][0] = 442.85614f;
    cgc.Value[587][1] = 0.0f;
    cgc.Value[588][0] = 442.86115f;
    cgc.Value[588][1] = 0.0f;
    cgc.Value[589][0] = 442.86615f;
    cgc.Value[589][1] = 0.0f;
    cgc.Value[590][0] = 442.87115f;
    cgc.Value[590][1] = 0.0f;
    cgc.Value[591][0] = 442.87613f;
    cgc.Value[591][1] = 0.0f;
    cgc.Value[592][0] = 442.88113f;
    cgc.Value[592][1] = 0.0f;
    cgc.Value[593][0] = 442.88614f;
    cgc.Value[593][1] = 0.0f;
    cgc.Value[594][0] = 442.89114f;
    cgc.Value[594][1] = 0.0f;
    cgc.Value[595][0] = 442.89615f;
    cgc.Value[595][1] = 0.0f;
    cgc.Value[596][0] = 442.90115f;
    cgc.Value[596][1] = 0.0f;
    cgc.Value[597][0] = 442.90613f;
    cgc.Value[597][1] = 0.0f;
    cgc.Value[598][0] = 442.91113f;
    cgc.Value[598][1] = 0.0f;
    cgc.Value[599][0] = 442.91614f;
    cgc.Value[599][1] = 0.0f;
    cgc.Value[600][0] = 442.92114f;
    cgc.Value[600][1] = 250000.0f;
    cgc.Value[601][0] = 442.92615f;
    cgc.Value[601][1] = 0.0f;
    cgc.Value[602][0] = 442.93115f;
    cgc.Value[602][1] = 0.0f;
    cgc.Value[603][0] = 442.93616f;
    cgc.Value[603][1] = 0.0f;
    cgc.Value[604][0] = 442.94113f;
    cgc.Value[604][1] = 0.0f;
    cgc.Value[605][0] = 442.94614f;
    cgc.Value[605][1] = 0.0f;
    cgc.Value[606][0] = 442.95114f;
    cgc.Value[606][1] = 0.0f;
    cgc.Value[607][0] = 442.95615f;
    cgc.Value[607][1] = 0.0f;
    cgc.Value[608][0] = 442.96115f;
    cgc.Value[608][1] = 0.0f;
    cgc.Value[609][0] = 442.96616f;
    cgc.Value[609][1] = 0.0f;
    cgc.Value[610][0] = 442.97113f;
    cgc.Value[610][1] = 0.0f;
    cgc.Value[611][0] = 442.97614f;
    cgc.Value[611][1] = 0.0f;
    cgc.Value[612][0] = 442.98114f;
    cgc.Value[612][1] = 0.0f;
    cgc.Value[613][0] = 442.98615f;
    cgc.Value[613][1] = 0.0f;
    cgc.Value[614][0] = 442.99115f;
    cgc.Value[614][1] = 0.0f;
    cgc.Value[615][0] = 442.99615f;
    cgc.Value[615][1] = 0.0f;
    cgc.Value[616][0] = 443.00113f;
    cgc.Value[616][1] = 0.0f;
    cgc.Value[617][0] = 443.00613f;
    cgc.Value[617][1] = 0.0f;
    cgc.Value[618][0] = 443.01114f;
    cgc.Value[618][1] = 0.0f;
    cgc.Value[619][0] = 443.01614f;
    cgc.Value[619][1] = 0.0f;
    cgc.Value[620][0] = 443.02115f;
    cgc.Value[620][1] = 0.0f;
    cgc.Value[621][0] = 443.02615f;
    cgc.Value[621][1] = 0.0f;
    cgc.Value[622][0] = 443.03113f;
    cgc.Value[622][1] = 0.0f;
    cgc.Value[623][0] = 443.03613f;
    cgc.Value[623][1] = 0.0f;
    cgc.Value[624][0] = 443.04114f;
    cgc.Value[624][1] = 83333.336f;
    cgc.Value[625][0] = 443.04614f;
    cgc.Value[625][1] = 0.0f;
    cgc.Value[626][0] = 443.05115f;
    cgc.Value[626][1] = 0.0f;
    cgc.Value[627][0] = 443.05615f;
    cgc.Value[627][1] = 0.0f;
    cgc.Value[628][0] = 443.06116f;
    cgc.Value[628][1] = 0.0f;
    cgc.Value[629][0] = 443.06613f;
    cgc.Value[629][1] = 0.0f;
    cgc.Value[630][0] = 443.07114f;
    cgc.Value[630][1] = 0.0f;
    cgc.Value[631][0] = 443.07614f;
    cgc.Value[631][1] = 0.0f;
    cgc.Value[632][0] = 443.08115f;
    cgc.Value[632][1] = 0.0f;
    cgc.Value[633][0] = 443.08615f;
    cgc.Value[633][1] = 0.0f;
    cgc.Value[634][0] = 443.09116f;
    cgc.Value[634][1] = 0.0f;
    cgc.Value[635][0] = 443.09613f;
    cgc.Value[635][1] = 0.0f;
    cgc.Value[636][0] = 443.10114f;
    cgc.Value[636][1] = 0.0f;
    cgc.Value[637][0] = 443.10614f;
    cgc.Value[637][1] = 0.0f;
    cgc.Value[638][0] = 443.11115f;
    cgc.Value[638][1] = 0.0f;
    cgc.Value[639][0] = 443.11615f;
    cgc.Value[639][1] = 0.0f;
    cgc.Value[640][0] = 443.12115f;
    cgc.Value[640][1] = 0.0f;
    cgc.Value[641][0] = 443.12613f;
    cgc.Value[641][1] = 0.0f;
    cgc.Value[642][0] = 443.13113f;
    cgc.Value[642][1] = 0.0f;
    cgc.Value[643][0] = 443.13614f;
    cgc.Value[643][1] = 0.0f;
    cgc.Value[644][0] = 443.14114f;
    cgc.Value[644][1] = 0.0f;
    cgc.Value[645][0] = 443.14615f;
    cgc.Value[645][1] = 0.0f;
    cgc.Value[646][0] = 443.15115f;
    cgc.Value[646][1] = 0.0f;
    cgc.Value[647][0] = 443.15613f;
    cgc.Value[647][1] = 0.0f;
    cgc.Value[648][0] = 443.16113f;
    cgc.Value[648][1] = 0.0f;
    cgc.Value[649][0] = 443.16614f;
    cgc.Value[649][1] = 0.0f;
    cgc.Value[650][0] = 443.17114f;
    cgc.Value[650][1] = 0.0f;
    cgc.Value[651][0] = 443.17615f;
    cgc.Value[651][1] = 0.0f;
    cgc.Value[652][0] = 443.18115f;
    cgc.Value[652][1] = 0.0f;
    cgc.Value[653][0] = 443.18616f;
    cgc.Value[653][1] = 0.0f;
    cgc.Value[654][0] = 443.19113f;
    cgc.Value[654][1] = 0.0f;
    cgc.Value[655][0] = 443.19614f;
    cgc.Value[655][1] = 0.0f;
    cgc.Value[656][0] = 443.20114f;
    cgc.Value[656][1] = 0.0f;
    cgc.Value[657][0] = 443.20615f;
    cgc.Value[657][1] = 0.0f;
    cgc.Value[658][0] = 443.21115f;
    cgc.Value[658][1] = 0.0f;
    cgc.Value[659][0] = 443.21616f;
    cgc.Value[659][1] = 0.0f;
    cgc.Value[660][0] = 443.22113f;
    cgc.Value[660][1] = 0.0f;
    cgc.Value[661][0] = 443.22614f;
    cgc.Value[661][1] = 0.0f;
    cgc.Value[662][0] = 443.23114f;
    cgc.Value[662][1] = 0.0f;
    cgc.Value[663][0] = 443.23615f;
    cgc.Value[663][1] = 0.0f;
    cgc.Value[664][0] = 443.24115f;
    cgc.Value[664][1] = 0.0f;
    cgc.Value[665][0] = 443.24615f;
    cgc.Value[665][1] = 0.0f;
    cgc.Value[666][0] = 443.25113f;
    cgc.Value[666][1] = 0.0f;
    cgc.Value[667][0] = 443.25613f;
    cgc.Value[667][1] = 0.0f;
    cgc.Value[668][0] = 443.26114f;
    cgc.Value[668][1] = 0.0f;
    cgc.Value[669][0] = 443.26614f;
    cgc.Value[669][1] = 0.0f;
    cgc.Value[670][0] = 443.27115f;
    cgc.Value[670][1] = 0.0f;
    cgc.Value[671][0] = 443.27615f;
    cgc.Value[671][1] = 0.0f;
    cgc.Value[672][0] = 443.28113f;
    cgc.Value[672][1] = 0.0f;
    cgc.Value[673][0] = 443.28613f;
    cgc.Value[673][1] = 0.0f;
    cgc.Value[674][0] = 443.29114f;
    cgc.Value[674][1] = 0.0f;
    cgc.Value[675][0] = 443.29614f;
    cgc.Value[675][1] = 0.0f;
    cgc.Value[676][0] = 443.30115f;
    cgc.Value[676][1] = 0.0f;
    cgc.Value[677][0] = 443.30615f;
    cgc.Value[677][1] = 0.0f;
    cgc.Value[678][0] = 443.31116f;
    cgc.Value[678][1] = 0.0f;
    cgc.Value[679][0] = 443.31613f;
    cgc.Value[679][1] = 0.0f;
    cgc.Value[680][0] = 443.32114f;
    cgc.Value[680][1] = 0.0f;
    cgc.Value[681][0] = 443.32614f;
    cgc.Value[681][1] = 0.0f;
    cgc.Value[682][0] = 443.33115f;
    cgc.Value[682][1] = 0.0f;
    cgc.Value[683][0] = 443.33615f;
    cgc.Value[683][1] = 0.0f;
    cgc.Value[684][0] = 443.34116f;
    cgc.Value[684][1] = 0.0f;
    cgc.Value[685][0] = 443.34613f;
    cgc.Value[685][1] = 0.0f;
    cgc.Value[686][0] = 443.35114f;
    cgc.Value[686][1] = 0.0f;
    cgc.Value[687][0] = 443.35614f;
    cgc.Value[687][1] = 0.0f;
    cgc.Value[688][0] = 443.36115f;
    cgc.Value[688][1] = 0.0f;
    cgc.Value[689][0] = 443.36615f;
    cgc.Value[689][1] = 0.0f;
    cgc.Value[690][0] = 443.37115f;
    cgc.Value[690][1] = 0.0f;
    cgc.Value[691][0] = 443.37613f;
    cgc.Value[691][1] = 0.0f;
    cgc.Value[692][0] = 443.38113f;
    cgc.Value[692][1] = 0.0f;
    cgc.Value[693][0] = 443.38614f;
    cgc.Value[693][1] = 0.0f;
    cgc.Value[694][0] = 443.39114f;
    cgc.Value[694][1] = 0.0f;
    cgc.Value[695][0] = 443.39615f;
    cgc.Value[695][1] = 0.0f;
    cgc.Value[696][0] = 443.40115f;
    cgc.Value[696][1] = 0.0f;
    cgc.Value[697][0] = 443.40613f;
    cgc.Value[697][1] = 0.0f;
    cgc.Value[698][0] = 443.41113f;
    cgc.Value[698][1] = 0.0f;
    cgc.Value[699][0] = 443.41614f;
    cgc.Value[699][1] = 0.0f;
    cgc.Value[700][0] = 443.42114f;
    cgc.Value[700][1] = 0.0f;
    cgc.Value[701][0] = 443.42615f;
    cgc.Value[701][1] = 0.0f;
    cgc.Value[702][0] = 443.43115f;
    cgc.Value[702][1] = 0.0f;
    cgc.Value[703][0] = 443.43616f;
    cgc.Value[703][1] = 0.0f;
    cgc.Value[704][0] = 443.44113f;
    cgc.Value[704][1] = 0.0f;
    cgc.Value[705][0] = 443.44614f;
    cgc.Value[705][1] = 0.0f;
    cgc.Value[706][0] = 443.45114f;
    cgc.Value[706][1] = 0.0f;
    cgc.Value[707][0] = 443.45615f;
    cgc.Value[707][1] = 0.0f;
    cgc.Value[708][0] = 443.46115f;
    cgc.Value[708][1] = 0.0f;
    cgc.Value[709][0] = 443.46616f;
    cgc.Value[709][1] = 0.0f;
    cgc.Value[710][0] = 443.47113f;
    cgc.Value[710][1] = 0.0f;
    cgc.Value[711][0] = 443.47614f;
    cgc.Value[711][1] = 0.0f;
    cgc.Value[712][0] = 443.48114f;
    cgc.Value[712][1] = 0.0f;
    cgc.Value[713][0] = 443.48615f;
    cgc.Value[713][1] = 0.0f;
    cgc.Value[714][0] = 443.49115f;
    cgc.Value[714][1] = 0.0f;
    cgc.Value[715][0] = 443.49615f;
    cgc.Value[715][1] = 0.0f;
    cgc.Value[716][0] = 443.50113f;
    cgc.Value[716][1] = 0.0f;
    cgc.Value[717][0] = 443.50613f;
    cgc.Value[717][1] = 0.0f;
    cgc.Value[718][0] = 443.51114f;
    cgc.Value[718][1] = 0.0f;
    cgc.Value[719][0] = 443.51614f;
    cgc.Value[719][1] = 0.0f;
    cgc.Value[720][0] = 443.52115f;
    cgc.Value[720][1] = 0.0f;
    cgc.Value[721][0] = 443.52615f;
    cgc.Value[721][1] = 0.0f;
    cgc.Value[722][0] = 443.53113f;
    cgc.Value[722][1] = 0.0f;
    cgc.Value[723][0] = 443.53613f;
    cgc.Value[723][1] = 0.0f;
    cgc.Value[724][0] = 443.54114f;
    cgc.Value[724][1] = 0.0f;
    cgc.Value[725][0] = 443.54614f;
    cgc.Value[725][1] = 0.0f;
    cgc.Value[726][0] = 443.55115f;
    cgc.Value[726][1] = 0.0f;
    cgc.Value[727][0] = 443.55615f;
    cgc.Value[727][1] = 0.0f;
    cgc.Value[728][0] = 443.56116f;
    cgc.Value[728][1] = 0.0f;
    cgc.Value[729][0] = 443.56613f;
    cgc.Value[729][1] = 0.0f;
    cgc.Value[730][0] = 443.57114f;
    cgc.Value[730][1] = 0.0f;
    cgc.Value[731][0] = 443.57614f;
    cgc.Value[731][1] = 0.0f;
    cgc.Value[732][0] = 443.58115f;
    cgc.Value[732][1] = 0.0f;
    cgc.Value[733][0] = 443.58615f;
    cgc.Value[733][1] = 0.0f;
    cgc.Value[734][0] = 443.59116f;
    cgc.Value[734][1] = 0.0f;
    cgc.Value[735][0] = 443.59613f;
    cgc.Value[735][1] = 0.0f;
    cgc.Value[736][0] = 443.60114f;
    cgc.Value[736][1] = 0.0f;
    cgc.Value[737][0] = 443.60614f;
    cgc.Value[737][1] = 0.0f;
    cgc.Value[738][0] = 443.61115f;
    cgc.Value[738][1] = 0.0f;
    cgc.Value[739][0] = 443.61615f;
    cgc.Value[739][1] = 0.0f;
    cgc.Value[740][0] = 443.62115f;
    cgc.Value[740][1] = 0.0f;
    cgc.Value[741][0] = 443.62613f;
    cgc.Value[741][1] = 0.0f;
    cgc.Value[742][0] = 443.63113f;
    cgc.Value[742][1] = 0.0f;
    cgc.Value[743][0] = 443.63614f;
    cgc.Value[743][1] = 0.0f;
    cgc.Value[744][0] = 443.64114f;
    cgc.Value[744][1] = 333333.34f;
    cgc.Value[745][0] = 443.64615f;
    cgc.Value[745][1] = 0.0f;
    cgc.Value[746][0] = 443.65115f;
    cgc.Value[746][1] = 0.0f;
    cgc.Value[747][0] = 443.65613f;
    cgc.Value[747][1] = 0.0f;
    cgc.Value[748][0] = 443.66113f;
    cgc.Value[748][1] = 0.0f;
    cgc.Value[749][0] = 443.66614f;
    cgc.Value[749][1] = 0.0f;
    cgc.Value[750][0] = 443.67114f;
    cgc.Value[750][1] = 0.0f;
    cgc.Value[751][0] = 443.67615f;
    cgc.Value[751][1] = 0.0f;
    cgc.Value[752][0] = 443.68115f;
    cgc.Value[752][1] = 0.0f;
    cgc.Value[753][0] = 443.68616f;
    cgc.Value[753][1] = 0.0f;
    cgc.Value[754][0] = 443.69113f;
    cgc.Value[754][1] = 0.0f;
    cgc.Value[755][0] = 443.69614f;
    cgc.Value[755][1] = 0.0f;
    cgc.Value[756][0] = 443.70114f;
    cgc.Value[756][1] = 0.0f;
    cgc.Value[757][0] = 443.70615f;
    cgc.Value[757][1] = 0.0f;
    cgc.Value[758][0] = 443.71115f;
    cgc.Value[758][1] = 0.0f;
    cgc.Value[759][0] = 443.71616f;
    cgc.Value[759][1] = 0.0f;
    cgc.Value[760][0] = 443.72113f;
    cgc.Value[760][1] = 0.0f;
    cgc.Value[761][0] = 443.72614f;
    cgc.Value[761][1] = 0.0f;
    cgc.Value[762][0] = 443.73114f;
    cgc.Value[762][1] = 0.0f;
    cgc.Value[763][0] = 443.73615f;
    cgc.Value[763][1] = 0.0f;
    cgc.Value[764][0] = 443.74115f;
    cgc.Value[764][1] = 0.0f;
    cgc.Value[765][0] = 443.74615f;
    cgc.Value[765][1] = 0.0f;
    cgc.Value[766][0] = 443.75113f;
    cgc.Value[766][1] = 0.0f;
    cgc.Value[767][0] = 443.75613f;
    cgc.Value[767][1] = 0.0f;
    cgc.Value[768][0] = 443.76114f;
    cgc.Value[768][1] = 0.0f;
    cgc.Value[769][0] = 443.76614f;
    cgc.Value[769][1] = 0.0f;
    cgc.Value[770][0] = 443.77115f;
    cgc.Value[770][1] = 0.0f;
    cgc.Value[771][0] = 443.77615f;
    cgc.Value[771][1] = 0.0f;
    cgc.Value[772][0] = 443.78113f;
    cgc.Value[772][1] = 0.0f;
    cgc.Value[773][0] = 443.78613f;
    cgc.Value[773][1] = 0.0f;
    cgc.Value[774][0] = 443.79114f;
    cgc.Value[774][1] = 0.0f;
    cgc.Value[775][0] = 443.79614f;
    cgc.Value[775][1] = 0.0f;
    cgc.Value[776][0] = 443.80115f;
    cgc.Value[776][1] = 0.0f;
    cgc.Value[777][0] = 443.80615f;
    cgc.Value[777][1] = 0.0f;
    cgc.Value[778][0] = 443.81116f;
    cgc.Value[778][1] = 0.0f;
    cgc.Value[779][0] = 443.81613f;
    cgc.Value[779][1] = 0.0f;
    cgc.Value[780][0] = 443.82114f;
    cgc.Value[780][1] = 0.0f;
    cgc.Value[781][0] = 443.82614f;
    cgc.Value[781][1] = 0.0f;
    cgc.Value[782][0] = 443.83115f;
    cgc.Value[782][1] = 0.0f;
    cgc.Value[783][0] = 443.83615f;
    cgc.Value[783][1] = 0.0f;
    cgc.Value[784][0] = 443.84116f;
    cgc.Value[784][1] = 0.0f;
    cgc.Value[785][0] = 443.84613f;
    cgc.Value[785][1] = 0.0f;
    cgc.Value[786][0] = 443.85114f;
    cgc.Value[786][1] = 0.0f;
    cgc.Value[787][0] = 443.85614f;
    cgc.Value[787][1] = 0.0f;
    cgc.Value[788][0] = 443.86115f;
    cgc.Value[788][1] = 0.0f;
    cgc.Value[789][0] = 443.86615f;
    cgc.Value[789][1] = 0.0f;
    cgc.Value[790][0] = 443.87115f;
    cgc.Value[790][1] = 0.0f;
    cgc.Value[791][0] = 443.87613f;
    cgc.Value[791][1] = 0.0f;
    cgc.Value[792][0] = 443.88113f;
    cgc.Value[792][1] = 166666.67f;
    cgc.Value[793][0] = 443.88614f;
    cgc.Value[793][1] = 0.0f;
    cgc.Value[794][0] = 443.89114f;
    cgc.Value[794][1] = 0.0f;
    cgc.Value[795][0] = 443.89615f;
    cgc.Value[795][1] = 0.0f;
    cgc.Value[796][0] = 443.90115f;
    cgc.Value[796][1] = 0.0f;
    cgc.Value[797][0] = 443.90613f;
    cgc.Value[797][1] = 0.0f;
    cgc.Value[798][0] = 443.91113f;
    cgc.Value[798][1] = 0.0f;
    cgc.Value[799][0] = 443.91614f;
    cgc.Value[799][1] = 0.0f;
    LipidomicsChromatogram cx = new LipidomicsChromatogram(cgc);
    cx.isProfile_ = true;
    cx.Mz = 10.973999f;
    cx.LowerMzBand = 3.973999f;
    cx.UpperMzBand = 17.973999f;
    
    SavGolJNI sav_gol_jni = new SavGolJNI();
    cx.Smooth(0f, 0, true, sav_gol_jni);
  }
  
  
  private void detectMSnWithLCB() {
    try {
//      LipidParameterSet param = new LipidParameterSet(464.546203613281f, "48", 1, "-OH", "", "C30 H59 N1 O3", "-H1 -O1",1,5);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 1993875328f;
//      probe1.AreaError = 398507968f;
//      probe1.Background = 4764792f;
//      probe1.Peak = 621.865234375f;
//      probe1.LowerValley = 611.065124511718f;
//      probe1.UpperValley = 631.46533203125f;
//      probe1.Mz = 464.351196289062f;
//      probe1.LowerMzBand = 0.43499755859375f;
//      probe1.UpperMzBand = 0.589996337890625f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
//      String[] chromPaths = StringUtils.getChromFilePaths("C:\\Sphingolipids\\dataStandardsPrelim\\20190111\\positive\\01092018 CER and SPH Mixes positive-10uM MIX 1.chrom");

//      LipidParameterSet param = new LipidParameterSet(829.798461914062f, "IS48", 1, "NH4", "", "C51 H89 O6 D7", "+N1 +H4",1,0);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 1993875328f;
//      probe1.AreaError = 398507968f;
//      probe1.Background = 4764792f;
//      probe1.Peak = 989.653991699218f;
//      probe1.LowerValley = 982.716003417968f;
//      probe1.UpperValley = 1004.40997314453f;
//      probe1.Mz = 829.798583984375f;
//      probe1.LowerMzBand = 0.0130615234375f;
//      probe1.UpperMzBand = 0.0159912109375f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
//      String[] chromPaths = StringUtils.getChromFilePaths("C:\\data\\Kristaps\\20171129\\TG quant NIST\\MCC007_Lipid01_NIST1_20171124_positive.chrom");

//      LipidParameterSet param = new LipidParameterSet(752.559936523437f, "38", 3, "-1", "", "C43 H80 N O7 P", "-H",1,0);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 1993875328f;
//      probe1.AreaError = 398507968f;
//      probe1.Background = 4764792f;
//      probe1.Peak = 1732.79003906258f;
//      probe1.LowerValley = 1706.17004394531f;
//      probe1.UpperValley = 1756.91003417968f;
//      probe1.Mz = 752.557983398437f;
//      probe1.LowerMzBand = 0.00897216796875f;
//      probe1.UpperMzBand = 0.00897216796875f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
//      String[] chromPaths = StringUtils.getChromFilePaths("C:\\data\\BiologicalExperiment\\Orbitrap_CID\\negative\\005_liver2-1_Orbitrap_CID_neg.chrom");

//      LipidParameterSet param = new LipidParameterSet(268.227111816406f, "16", 1, "-OH", "", "C16 H31 N O3", "-H1 -O1",1,2);
//
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 1993875328f;
//      probe1.AreaError = 398507968f;
//      probe1.Background = 4764792f;
//      probe1.Peak = 82.083999633789f;
//      probe1.LowerValley = 63.7661018371582f;
//      probe1.UpperValley = 105.580001831054f;
//      probe1.Mz = 268.228088378906f;
//      probe1.LowerMzBand = 0.003997802734375f;
//      probe1.UpperMzBand = 0.0050048828125f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
//      String[] chromPaths = StringUtils.getChromFilePaths("C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\Mix1\\Mix1_1.chrom");

//      LipidParameterSet param = new LipidParameterSet(442.090179443359f, "26", 0, "-H", "", "C26 H53 N1 O4", "-H",1,3);
//      CgProbe probe1 = new CgProbe(0,1);
//      probe1.AreaStatus = CgAreaStatus.OK;
//      probe1.Area = 199908520f;
//      probe1.AreaError = 398507968f;
//      probe1.Background = 4764792f;
//      probe1.Peak = 504.243072509765f;
//      probe1.LowerValley = 494.642974853515f;
//      probe1.UpperValley = 513.843139648437f;
//      probe1.Mz = 442.265167236328f;
//      probe1.LowerMzBand = 0.504974365234375f;
//      probe1.UpperMzBand = 0.58502197265625f;
//      probe1.isotopeNumber = 0;
//      param.AddProbe(probe1);
//      String[] chromPaths = StringUtils.getChromFilePaths("C:\\Sphingolipids\\dataStandardsPrelim\\20190205\\negative\\02042019 IDA neg new-new CER SPH Mix 1 neg.chrom");

      LipidParameterSet param = new LipidParameterSet(600.520874023437f, "34", 0, "HCOO", "", "C34 H69 N O4", "C1 H1 O2",1,3);

      CgProbe probe1 = new CgProbe(0,1);
      probe1.AreaStatus = CgAreaStatus.OK;
      probe1.Area = 1993875328f;
      probe1.AreaError = 398507968f;
      probe1.Background = 4764792f;
      probe1.Peak = 1244.77001953125f;
      probe1.LowerValley = 1227.31994628906f;
      probe1.UpperValley = 1259.92004394531f;
      probe1.Mz = 600.52294921875f;
      probe1.LowerMzBand = 0.00506591796875f;
      probe1.UpperMzBand = 0.00897216796875f;
      probe1.isotopeNumber = 0;
      param.AddProbe(probe1);
      String[] chromPaths = StringUtils.getChromFilePaths("C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\Mix1\\Mix1_neg_1.chrom");

      
      System.out.println(chromPaths[1]);
      LipidomicsAnalyzer lAnalyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
      setStandardParameters(lAnalyzer);
      MSnAnalyzer analyzer = new MSnAnalyzer(null,"Cer","HCOO",param,lAnalyzer,null,false,true,true);
//      MSnAnalyzer analyzer = new MSnAnalyzer(null,"P-PE","-H",param,lAnalyzer,null,false,true,true);
//      MSnAnalyzer analyzer = new MSnAnalyzer(null,"TG","NH4",param,lAnalyzer,null,false,true,true);
      
      System.out.println("Status: "+analyzer.checkStatus());
      if (!(analyzer.getResult() instanceof LipidomicsMSnSet)) return;
      LipidomicsMSnSet result = (LipidomicsMSnSet)analyzer.getResult();
      for (Object nameObject : result.getMSnIdentificationNames()){
        if (nameObject instanceof Vector){
          for (String name : (Vector<String>)nameObject){
            System.out.println(result.getRt()+": "+name);
          }
        }else{
          String name = (String)nameObject;
          System.out.println(result.getNameIncludingModification()+": "+name);
        }
      }

      if (analyzer.checkStatus()!=LipidomicsMSnSet.NO_MSN_PRESENT && analyzer.checkStatus()!=LipidomicsMSnSet.DISCARD_HIT){
//        Hashtable<String,CgProbe> headFragments = analyzer.getHeadGroupFragments();
//        for (String fragmentName : headFragments.keySet()){
//          CgProbe probe = headFragments.get(fragmentName);
//          System.out.println(fragmentName+": \t"+probe.Mz+"\t"+probe.Area);
//        }
//        for (IntensityRuleVO fulfilled : analyzer.getFulfilledHeadIntensityRules().values()) {
//          System.out.println("fulfilled: "+fulfilled.getRuleIdentifier());
//        }
        Hashtable<String,Hashtable<String,CgProbe>> chainFragments = analyzer.getChainFragments();
        for (String chain : chainFragments.keySet()){
          System.out.println("1. "+chain);
          Hashtable<String,CgProbe> frags = chainFragments.get(chain);
          for (String frag : frags.keySet()){
            CgProbe probe = frags.get(frag);
            System.out.println(frag+": \t"+probe.Mz+"\t"+probe.Area);
          }
        }
        for (String key : analyzer.getFulfilledChainIntensityRules().keySet()) {
          Hashtable<String,IntensityChainVO> rules = analyzer.getFulfilledChainIntensityRules().get(key);
          for (IntensityChainVO rule : rules.values()) {
            System.out.println("fulfilled: "+key+"   "+rule.getRuleIdentifier());
          }
        }
      }
    }catch(Exception ex) {
      ex.printStackTrace();
    }
  }
  
  private void parseMSDialTxt() {
    //MSDialTxtParser dialParser = new MSDialTxtParser("C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\MS-Dial_Mix1\\export\\Mix1_neg_1.txt");
    MSDialTxtParser dialParser = new MSDialTxtParser("C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\MS-Dial_Mix1\\export\\Mix1_1.txt");
    try {
      dialParser.parse();
      Vector<MSDialEntry> entries = dialParser.getResults();
//      for (MSDialEntry entry : entries) {
//        System.out.println(entry.getName());
//      }
      
    }
    catch (Exception e) {
      e.printStackTrace();
    }
    
  }
  
  
  private void compareLDABMSDialControlledPositiveProbes() {
    LinkedHashMap<String,LinkedHashMap<String,Boolean[]>> comparableClassesAndAdducts = new LinkedHashMap<String,LinkedHashMap<String,Boolean[]>>();
    Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>> correctAnalytes = new Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>>();
    Hashtable<String,Integer> snPositions = new Hashtable<String,Integer>();
    
    String lipidClass = "SphBase";
    //the first boolean tells if this adduct should be included in statistics, the second if it the analytes should be RT-filtered, the third one if chain identification should be included in statistics
    LinkedHashMap<String,Boolean[]> adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true,true});
    adducts.put("Na",new Boolean[]{false,true,true});
    adducts.put("-OH",new Boolean[]{false,true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 1);
    //correctAnalytes.put(lipidClass, Mix1Standards.getSphBaseStandards());
    correctAnalytes.put(lipidClass, Mix2Standards.getSphBaseStandards());
    
    lipidClass = "Cer";
    //the first boolean tells if this adduct should be included in statistics, the second if it the analytes should be RT-filtered, the third one if chain identification should be included in statistics
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true,true});
    adducts.put("Na",new Boolean[]{true,true,false});
    adducts.put("-OH",new Boolean[]{false,true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    //correctAnalytes.put(lipidClass, Mix1Standards.getCerStandards());
    correctAnalytes.put(lipidClass, Mix2Standards.getCerStandards());
    
    lipidClass = "Cer1P";
    //the first boolean tells if this adduct should be included in statistics, the second if it the analytes should be RT-filtered, the third one if chain identification should be included in statistics
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,false,true});
    adducts.put("Na",new Boolean[]{true,false,false});
    adducts.put("-OH",new Boolean[]{false,false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    //correctAnalytes.put(lipidClass, Mix1Standards.getCer1PStandards());
    correctAnalytes.put(lipidClass, Mix2Standards.getCer1PStandards());
    
    lipidClass = "HexCer";
    //the first boolean tells if this adduct should be included in statistics, the second if it the analytes should be RT-filtered, the third one if chain identification should be included in statistics
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true,true});
    adducts.put("Na",new Boolean[]{true,true,false});
    adducts.put("-OH",new Boolean[]{false,true,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    //correctAnalytes.put(lipidClass, Mix1Standards.getHexCerStandards());
    correctAnalytes.put(lipidClass, Mix2Standards.getHexCerStandards());
    
    lipidClass = "SM";
    //the first boolean tells if this adduct should be included in statistics, and the second if it the analytes should be RT-filtered
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("H",new Boolean[]{true,true,false});
    adducts.put("Na",new Boolean[]{true,true,false});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    //correctAnalytes.put(lipidClass, Mix1Standards.getSMStandards());
    correctAnalytes.put(lipidClass, Mix2Standards.getSMStandards());

    //correcting the RT shift in positive ion mode
    for (LinkedHashMap<String,ReferenceInfoVO> classCorrect : correctAnalytes.values()) {
      for (ReferenceInfoVO info : classCorrect.values()) {
        for (int i=0; i!=info.getCorrectRts().length; i++) {
          info.getCorrectRts()[i] = info.getCorrectRts()[i]-0.15d;        }
      }
    }
    
//    String chromFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\Mix1\\Mix1_5.chrom";
//    String ldaFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\Mix1\\Mix1_5_stands_positive_Mix1.xlsx";
//    String msDialFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\MS-Dial_Mix1\\export\\Mix1_5.txt";
//    String quantFile = "C:\\Sphingolipids\\Experiment1\\massLists\\Orbitrap\\stands_positive_Mix1.xlsx";
//    String outputFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\MS-Dial_Mix1\\Mix1_5_Dial_comp_generated.xlsx";

    String chromFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\Mix2\\Mix2_5.chrom";
    String ldaFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\Mix2\\Mix2_5_stands_positive_Mix2.xlsx";
    String msDialFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\MS-Dial_Mix2\\export\\Mix2_5.txt";
    String quantFile = "C:\\Sphingolipids\\Experiment1\\massLists\\Orbitrap\\stands_positive_Mix2.xlsx";
    String outputFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\positive\\MS-Dial_Mix2\\Mix2_5_Dial_comp_generated.xlsx";

    
    compareLDAMSDialControlledProbes(comparableClassesAndAdducts,snPositions,correctAnalytes,chromFile,quantFile,ldaFile,msDialFile,/*msFinderFile,*/outputFile);

    
  }
  
  
  private void compareLDABMSDialControlledNegativeProbes() {
    LinkedHashMap<String,LinkedHashMap<String,Boolean[]>> comparableClassesAndAdducts = new LinkedHashMap<String,LinkedHashMap<String,Boolean[]>>();
    Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>> correctAnalytes = new Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>>();
    Hashtable<String,Integer> snPositions = new Hashtable<String,Integer>();
    
    String lipidClass = "Cer";
    //the first boolean tells if this adduct should be included in statistics, the second if it the analytes should be RT-filtered, the third one if chain identification should be included in statistics
    LinkedHashMap<String,Boolean[]> adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("HCOO",new Boolean[]{true,true,false});
    adducts.put("-H",new Boolean[]{true,true,true});
    adducts.put("Cl",new Boolean[]{false,true,false});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    //correctAnalytes.put(lipidClass, Mix1Standards.getCerStandards());
    correctAnalytes.put(lipidClass, Mix2Standards.getCerStandards());
    
    lipidClass = "Cer1P";
    //the first boolean tells if this adduct should be included in statistics, the second if it the analytes should be RT-filtered, the third one if chain identification should be included in statistics
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("-H",new Boolean[]{true,false,true});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    //correctAnalytes.put(lipidClass, Mix1Standards.getCer1PStandards());
    correctAnalytes.put(lipidClass, Mix2Standards.getCer1PStandards());
    
    lipidClass = "HexCer";
    //the first boolean tells if this adduct should be included in statistics, the second if it the analytes should be RT-filtered, the third one if chain identification should be included in statistics
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("HCOO",new Boolean[]{true,true,false});
    adducts.put("-H",new Boolean[]{true,true,true});
    adducts.put("Cl",new Boolean[]{false,true,false});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    //correctAnalytes.put(lipidClass, Mix1Standards.getHexCerStandards());
    correctAnalytes.put(lipidClass, Mix2Standards.getHexCerStandards());
    
    lipidClass = "SM";
    //the first boolean tells if this adduct should be included in statistics, and the second if it the analytes should be RT-filtered
    adducts = new LinkedHashMap<String,Boolean[]>();
    adducts.put("HCOO",new Boolean[]{true,true,false});
    adducts.put("Cl",new Boolean[]{false,true,false});
    comparableClassesAndAdducts.put(lipidClass, adducts);
    snPositions.put(lipidClass, 2);
    //correctAnalytes.put(lipidClass, Mix1Standards.getSMStandards());
    correctAnalytes.put(lipidClass, Mix2Standards.getSMStandards());

    
//    String chromFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\Mix1\\Mix1_neg_5.chrom";
//    String ldaFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\Mix1\\Mix1_neg_5_stands_negative_Mix1.xlsx";
//    String msDialFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\MS-Dial_Mix1\\export\\Mix1_neg_5.txt";
//    String quantFile = "C:\\Sphingolipids\\Experiment1\\massLists\\Orbitrap\\stands_negative_Mix1.xlsx";
    //String msFinderFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\MS-Dial_Mix1\\export_MSFinder2\\Structure result-2101.txt";
    //String msFinderFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\MS-Dial_Mix1\\export_MSFinder\\Structure result-2097.txt";
//    String outputFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\MS-Dial_Mix1\\Mix1_neg_5_Dial_comp_generated.xlsx";

    String chromFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\Mix2\\Mix2_neg_5.chrom";
    String ldaFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\Mix2\\Mix2_neg_5_stands_negative_Mix2.xlsx";
    String msDialFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\MS-Dial_Mix2\\export\\Mix2_neg_5.txt";
    String quantFile = "C:\\Sphingolipids\\Experiment1\\massLists\\Orbitrap\\stands_negative_Mix2.xlsx";
    String outputFile = "C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\MS-Dial_Mix2\\Mix2_neg_5_Dial_comp_generated.xlsx";
    
    compareLDAMSDialControlledProbes(comparableClassesAndAdducts,snPositions,correctAnalytes,chromFile,quantFile,ldaFile,msDialFile,/*msFinderFile,*/outputFile);

  }
  
  private void compareLDAMSDialControlledProbes(LinkedHashMap<String,LinkedHashMap<String,Boolean[]>> comparableClassesAndAdducts,
      Hashtable<String,Integer> snPositions, Hashtable<String,LinkedHashMap<String,ReferenceInfoVO>> correctAnalytes,
      String chromFile, String quantFile, String ldaFile, String msdialFile,/* String msfinderFile,*/ String outputFile){
    try{
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outputFile));
      Workbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      CellStyle leftHeaderStyle = resultWorkbook.createCellStyle();
      leftHeaderStyle.cloneStyleFrom(headerStyle);
      leftHeaderStyle.setAlignment(CellStyle.ALIGN_LEFT);
      headerStyle = leftHeaderStyle;
      CellStyle centerStyle = getCenterStyle(resultWorkbook);

      CellStyle fullyCorrectStyle = getFullyCorrectStyle(resultWorkbook);
      CellStyle faCorrectStyle = getFACorrectStyle(resultWorkbook);
      CellStyle faFoundStyle = getFAFoundStyle(resultWorkbook);
      CellStyle ms1FoundStyle = getMS1FoundStyle(resultWorkbook);
      CellStyle notFoundStyle = getNotFoundStyle(resultWorkbook);
      Hashtable<String,Vector<LipidParameterSet>> resultsLDA = LDAResultReader.readResultFile(ldaFile, new Hashtable<String,Boolean>()).getIdentifications();
      MSDialTxtParser msdialParser = new MSDialTxtParser(msdialFile);
      msdialParser.parse();
      //MSFinderStructureParser finderParser = new MSFinderStructureParser(msfinderFile);
      //finderParser.parse();
      //TODO: this might need to be changed
     // Hashtable<String,Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>> resultsLB = msdialParser.getResults();
      Vector<MSDialEntry> resultsDial = msdialParser.getResults();
      //Vector<MSFinderEntry> resultsFinder = finderParser.getResults();
      String[] chromPaths = StringUtils.getChromFilePaths(chromFile);   
      LipidomicsAnalyzer analyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantVOs = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>)QuantificationThread.parseQuantExcelFile(quantFile,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f,true).get(3);

      ////Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantVOs = null;
      Sheet summarySheet = resultWorkbook.createSheet("Summary");
      Sheet detailsSheet = resultWorkbook.createSheet("Details");

      int ms1Detectables = 0;
      int faDetectables = 0;
      int positionDetectables = 0;
      
      int ms1DetectedLDA = 0;
      int faDetectedLDAWOFPs = 0;
      int faDetectedLDA = 0;

      int ms1DetectedDial = 0;
      int faDetectedDialWOFPs = 0;
      int faDetectedDial = 0;
      
      int rowCountDetails = 0;
      int rowCountSummary = 0;
      
      int classColumn = 0;
      int correctColumn = 1;
      int width = (12+1)*256;
      detailsSheet.setColumnWidth(correctColumn, width);
      int adductColumn = 2;
      int ldaColumn = 3;
      width = (26+1)*256;
      detailsSheet.setColumnWidth(ldaColumn, width);
      int dialColumn = 4;
      detailsSheet.setColumnWidth(dialColumn, width);
      int dialProb = 5;
      int ldaRt = 6;
/****      int lbRts = 7;
      int commentColumn = 8;*/
      int commentColumn = 7;
      
      int sumTypeColumn = 0;
      int ldaMs1Column = 1;
      int ldaMs1PercentColumn = 2;
      int dialMs1Column = 3;
      int dialMs1PercentColumn = 4;
      int ldaFaColumn = 5;
      int ldaFaPercentColumn = 6;
      int dialFaColumn = 7;
      int dialFaPercentColumn = 8;
      int ldaFaWoFPColumn = 9;
      int ldaFaWoFPPercentColumn = 10;
      int dialFaWoFPColumn = 11;
      int dialFaWoFPPercentColumn = 12;

      Row sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      Cell sumCell = sumRow.createCell(ldaMs1Column);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA identified");
      sumCell = sumRow.createCell(dialMs1Column);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("MS-DIAL identified");
      sumCell = sumRow.createCell(ldaFaColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA FAs");
      sumCell = sumRow.createCell(dialFaColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("MS-DIAL FAs");
      sumCell = sumRow.createCell(ldaFaWoFPColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA FAs no FPs");
      sumCell = sumRow.createCell(dialFaWoFPColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("MS-DIAL FAs no FPs");
      int typeColumn = 0;
      int ldaHitsColumn = 4;
      int ldaPercentColumn = 5;
      int lbHitsColumn = 6;
      int lbPercentColumn = 7;

      
      LinkedHashMap<String,LinkedHashMap<String,LdaDialStandardsEvidence>> speciesEvidence = new LinkedHashMap<String,LinkedHashMap<String,LdaDialStandardsEvidence>>();
      
      for (String lipidClass : comparableClassesAndAdducts.keySet()){
        LinkedHashMap<String,Boolean[]> comparableAdducts = comparableClassesAndAdducts.get(lipidClass);
        int numberSns = snPositions.get(lipidClass);
        LinkedHashMap<String,ReferenceInfoVO> ms1Analytes = correctAnalytes.get(lipidClass);
        Vector<LipidParameterSet> ldaAnalytes = resultsLDA.get(lipidClass);
        //TODO: this part is for testing:
        System.out.println("!!!!!!!!!!!!!!!!!!! "+lipidClass+" !!!!!!!!!!!!!!!!!!!");
        /****        for (String name : ms1Analytes.keySet()) {
          ReferenceInfoVO info = ms1Analytes.get(name);
          for (String adduct : info.getAdducts().keySet()) {
            float mz = (float)info.getAdducts().get(adduct).doubleValue();
            String rts = "";
            for (double rt : info.getCorrectRts()) {
              if (rts.length()>0)
                rts += ";";
              rts += Calculator.FormatNumber(rt, 2);
            }
            System.out.println("1. "+name+"_"+adduct+"_"+rts);
            for (MSDialEntry entry : resultsDial) {
              if ((mz-0.006d)<entry.getMz() && entry.getMz()<(mz+0.006d)){
                boolean found = false;
                for (int i=0; i!=info.getCorrectRts().length; i++) {
                  float rt = (float)info.getCorrectRts()[i];
                  if (lipidClass.equalsIgnoreCase("Cer1P") || ((rt-0.25f)<entry.getRt() && entry.getRt()<(rt+0.25f))) {
                    found = true;
                  }
                }
                if (found)
                  System.out.println("MS-Dial: "+entry.getName()+" "+entry.getRt()+" "+entry.getScore());
              }
            }
            for (MSFinderEntry entry : resultsFinder) {
              if ((mz-0.006d)<entry.getMz() && entry.getMz()<(mz+0.006d)){
                boolean found = false;
                for (int i=0; i!=info.getCorrectRts().length; i++) {
                  float rt = (float)info.getCorrectRts()[i];
                  if (lipidClass.equalsIgnoreCase("Cer1P") || ((rt-0.25f)<entry.getRt() && entry.getRt()<(rt+0.25f))) {
                    found = true;
                  }
                }
                if (found) {
                  Hashtable<Integer,MSFinderHitVO> hits = entry.getHits();
                  for (int i=0; i!=hits.size(); i++) {
                    MSFinderHitVO hit = hits.get((i+1));
                    System.out.println("MS-Finder: "+(i+1)+" "+hit.getStructureName()+"_"+entry.getAdduct()+" ; "+entry.getRt()+" ; "+hit.getScore());
                    
                  }                  
                }

              }
            }
          }*/
        ////}
              
        LinkedHashMap<String,LdaDialStandardsEvidence> speciesEvidenceClass = new LinkedHashMap<String,LdaDialStandardsEvidence>();
        Row row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;
        Cell cell = row.createCell(classColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Class");
        cell = row.createCell(correctColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Correct");
        cell = row.createCell(adductColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Adduct");
        cell = row.createCell(ldaColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA");
        cell = row.createCell(dialColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("MS-DIAL");
        cell = row.createCell(dialProb);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("DIAL-Prob");
        cell = row.createCell(ldaRt);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA-RT");
/****        cell = row.createCell(lbRts);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LB-RTs");*/
        cell = row.createCell(commentColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Comments");

        int ms1DetectablesClass = 0;
        int faDetectablesClass = 0;
        int positionDetectablesClass = 0;
        
        int ms1DetectedLDAClass = 0;
        int faDetectedLDAWOFPsClass = 0;
        int faDetectedLDAClass = 0;

        int ms1DetectedDialClass = 0;
        int faDetectedDialWOFPsClass = 0;
        int faDetectedDialClass = 0;
        
        for (String adduct : comparableAdducts.keySet()){

          int ms1DetectablesAdduct = 0;
          int faDetectablesAdduct = 0;
          int positionDetectablesAdduct = 0;
          
          int ms1DetectedLDAAdduct = 0;
          int faDetectedLDAWOFPsAdduct = 0;
          int faDetectedLDAAdduct = 0;

          int ms1DetectedDialAdduct = 0;
          int faDetectedDialWOFPsAdduct = 0;
          int faDetectedDialAdduct = 0;

          boolean includeInAnalysis = comparableAdducts.get(adduct)[0];
          boolean includeChains = comparableAdducts.get(adduct)[2];
          String commentText;
          
          for (String analyte : ms1Analytes.keySet()){
            ReferenceInfoVO correctStructureInfo = ms1Analytes.get(analyte);
            String rts = "";
            commentText = "";
            if (!includeChains) 
              commentText = "chains excluded from statistics";
            for (double rt : correctStructureInfo.getCorrectRts()) {
              if (rts.length()>0)
                rts += ";";
              rts += Calculator.FormatNumber(rt, 2);
            }

            //System.out.println("1. "+correctStructureInfo.getMS2Name()+"_"+adduct+"_"+rts);
            LipidomicsMSnSet ldaAnalyte = (LipidomicsMSnSet)getLDAAnalyte(ldaAnalytes,analyte,adduct,false,lipidClass.equalsIgnoreCase("Cer1P"),correctStructureInfo.getCorrectRts());
//            if (correctStructureInfo==null)
//              System.out.println("11111");
//            if (correctStructureInfo.getAdducts()==null)
//              System.out.println("22222");
//            if (correctStructureInfo.getAdducts().get(adduct)==null)
//              System.out.println("33333");

            Vector<Vector<MSDialEntry>> inAndoutRt = getDialAnalytes(resultsDial, correctStructureInfo.getAdducts().get(adduct).doubleValue(), correctStructureInfo.getCorrectRts(),comparableAdducts.get(adduct)[1]);
            Vector<MSDialEntry> filteredDial = inAndoutRt.get(0);
            Vector<MSDialEntry> outRt = inAndoutRt.get(1);
            //            for (MSDialEntry entry : filteredDial) {
//              System.out.println("MS-Dial: "+entry.getName()+" "+entry.getRt()+" "+entry.getScore());
//            }
            //LipidBLASTIdentificationVO identVO = getDialAnalyte(lbAnalytes,analyte,adduct,correctStructureInfo.getCorrectRts()[0],comparableAdducts.get(adduct)[1]);
            String correctStructure = correctStructureInfo.getMS2Name();
            
            LdaDialStandardsEvidence ldaDialEvidence = new LdaDialStandardsEvidence(analyte);
            if (speciesEvidenceClass.containsKey(analyte)) ldaDialEvidence = speciesEvidenceClass.get(analyte);
              row = detailsSheet.createRow(rowCountDetails);
              rowCountDetails++;
              cell = row.createCell(classColumn);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(lipidClass);
              cell = row.createCell(correctColumn);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(correctStructure);
              cell = row.createCell(adductColumn);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(adduct);
              
              Object[] lda = checkLDAEvidence(correctStructure, ldaAnalyte);
              String nameString = (String)lda[1];
              int evidence = (Integer)lda[0];
              //System.out.println(nameString+"_"+adduct+": "+evidence);
                            
              Object[] msDial = checkDialEvidence(lipidClass, analyte, correctStructure, filteredDial, numberSns);
              String dialString = (String)msDial[1];
              int dialEvidence = (Integer)msDial[0];
              //System.out.println(nameString+"_"+adduct+": "+dialEvidence);
              
              boolean msMSPresent = true;
              LipidParameterSet ldaMS1Result = null;
              if (evidence==0 && dialEvidence==0){
                ldaMS1Result = getLDAAnalyte(ldaAnalytes,analyte,adduct,true,true,correctStructureInfo.getCorrectRts());
                if (!areThereAnyMSnSpectra(analyzer,ldaMS1Result,quantVOs.get(lipidClass).get(analyte).get(adduct))){
                  evidence = -1;
                  nameString = "no MS/MS";
                  if (ldaMS1Result!=null) nameString = ldaMS1Result.getNameStringWithoutRt();      
                  dialEvidence = -1;
                  dialString = "no MS/MS";
                  msMSPresent = false;
                }               
              }
              cell = row.createCell(ldaColumn);
              CellStyle style = getStyleBasedOnEvidenceDial(evidence, fullyCorrectStyle, faCorrectStyle, faFoundStyle,
                  ms1FoundStyle, notFoundStyle);
              if (evidence==0 && dialEvidence>0 ){
                ldaMS1Result = getLDAAnalyte(ldaAnalytes,analyte,adduct,true,true,correctStructureInfo.getCorrectRts());
                if (ldaMS1Result!=null){
                  nameString = ldaMS1Result.getNameStringWithoutRt();
                  style=null;
                }
              }
              if (style!=null) cell.setCellStyle(style);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(nameString);

              cell = row.createCell(dialColumn);
              style = getStyleBasedOnEvidenceDial(dialEvidence, fullyCorrectStyle, faCorrectStyle, faFoundStyle,
                  ms1FoundStyle, notFoundStyle);
              if (style!=null) cell.setCellStyle(style);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              cell.setCellValue(dialString);

              if (dialEvidence>0){
                cell = row.createCell(dialProb);
                cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                cell.setCellValue((Double)msDial[2]);
                
              }
              if (evidence>0 || ldaMS1Result!=null){
                cell = row.createCell(ldaRt);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                if (evidence>0)
                  cell.setCellValue(ldaAnalyte.getRt());
                else
                  cell.setCellValue(ldaMS1Result.getRt());
              }
              if (includeInAnalysis){
                if (msMSPresent && correctStructureInfo.useInEvaluation()){
                  ms1DetectablesClass++;
                  ms1DetectablesAdduct++;
                  //for debugging
//                  if ((blastEvidence>=FA_FOUND)){
//                    cell = row.createCell(commentColumn);
//                    cell.setCellType(Cell.CELL_TYPE_STRING);
//                    cell.setCellValue(((Vector<String>)lBlast[3]).toString());                
//                  }
                  if (dialEvidence==0 && outRt.size()>0){
                    //WAS-DONE: for backward compatibility, the following was changed: getCorrectRts()[0]
                    String startRt = Calculator.FormatNumberToString(correctStructureInfo.getCorrectRts()[0]-MSDIAL_RT_TOL, 2);
                    String stopRt = Calculator.FormatNumberToString(correctStructureInfo.getCorrectRts()[0]+MSDIAL_RT_TOL, 2);
                    String outRts = "";
                    for (MSDialEntry outHit : outRt) outRts+=outHit.getName()+" "+Calculator.FormatNumber(outHit.getRt(),2d)+";";
                    outRts = outRts.substring(0,outRts.length()-1);
                    if (commentText.length()>0)
                      commentText += "; ";
                    commentText += "MS-DIAL reported hits outside the RT-tolerance "+startRt+"-"+stopRt+"min: "+outRts;
                    
                  }
                  if (commentText.length()>0) {
                    cell = row.createCell(commentColumn);
                    cell.setCellType(Cell.CELL_TYPE_STRING);
                    cell.setCellValue(commentText);
                  }
                  if (evidence>=MS1_FOUND){
                    ms1DetectedLDAClass++;
                    ms1DetectedLDAAdduct++;
                    ldaDialEvidence.setMs1DetectedLDA(true);
                  }
                  if (dialEvidence>=MS1_FOUND){
                    ms1DetectedDialClass++;
                    ms1DetectedDialAdduct++;
                    ldaDialEvidence.setMs1DetectedLB(true);
                  }
                  if (numberSns>1 && includeChains){
                    faDetectablesClass++;
                    faDetectablesAdduct++;
                    ldaDialEvidence.setFaDetectable(true);
                    if (correctStructureInfo.isPositionAvailable() && !(lipidClass.equalsIgnoreCase("DG") && adduct.equalsIgnoreCase("NH4"))){
                      positionDetectablesClass++;
                      positionDetectablesAdduct++;
                    }
                    if (evidence>=FA_CORRECT){
                      faDetectedLDAWOFPsClass++;
                      faDetectedLDAWOFPsAdduct++;
                      ldaDialEvidence.setFaDetectedLDAWoFPs(true);
                    }
                    if (dialEvidence>=FA_CORRECT){
                      faDetectedDialWOFPsClass++;
                      faDetectedDialWOFPsAdduct++;
                      ldaDialEvidence.setFaDetectedLBWoFPs(true);
                    }
                    if (evidence>=FA_FOUND){
                      faDetectedLDAClass++;
                      faDetectedLDAAdduct++;
                      ldaDialEvidence.setFaDetectedLDA(true);
                    }
                    if (dialEvidence>=FA_FOUND){
                      faDetectedDialClass++;
                      faDetectedDialAdduct++;
                      ldaDialEvidence.setFaDetectedLB(true);
                    }
                  }
                  speciesEvidenceClass.put(analyte,ldaDialEvidence);
//                } else if (!correctStructures.get(correctStructure)[0]){
                }else{  
                  cell = row.createCell(commentColumn);
                  cell.setCellType(Cell.CELL_TYPE_STRING);
                  cell.setCellValue("excluded from statistics");                  
                }
              }else{
                cell = row.createCell(commentColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue("excluded from statistics");                
              }
//            }
          }
          
          if (ms1DetectablesAdduct==0) continue;
          
          sumRow = summarySheet.createRow(rowCountSummary);
          rowCountSummary++;
          sumCell = sumRow.createCell(sumTypeColumn);
          sumCell.setCellStyle(leftHeaderStyle);
          sumCell.setCellValue(lipidClass+"_"+adduct);
          
          sumCell = sumRow.createCell(ldaMs1Column);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          sumCell.setCellValue(ms1DetectedLDAAdduct+"/"+ms1DetectablesAdduct);
          sumCell = sumRow.createCell(ldaMs1PercentColumn);
          sumCell.getCellStyle().setAlignment(CellStyle.ALIGN_CENTER);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*ms1DetectedLDAAdduct))/((double)ms1DetectablesAdduct), 0)+"%");
          sumCell = sumRow.createCell(dialMs1Column);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          sumCell.setCellValue(ms1DetectedDialAdduct+"/"+ms1DetectablesAdduct);
          sumCell = sumRow.createCell(dialMs1PercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*ms1DetectedDialAdduct))/((double)ms1DetectablesAdduct), 0)+"%");
          sumCell = sumRow.createCell(ldaFaColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (includeChains && numberSns>1) sumCell.setCellValue(faDetectedLDAAdduct+"/"+faDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(ldaFaPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (includeChains && numberSns>1 && faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLDAAdduct))/((double)faDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(dialFaColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (includeChains && numberSns>1 && faDetectablesClass>0) sumCell.setCellValue(faDetectedDialAdduct+"/"+faDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(dialFaPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (includeChains && numberSns>1 && faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedDialAdduct))/((double)faDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(ldaFaWoFPColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (includeChains && numberSns>1 && faDetectablesClass>0) sumCell.setCellValue(faDetectedLDAWOFPsAdduct+"/"+faDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(ldaFaWoFPPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (includeChains && numberSns>1 && faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLDAWOFPsAdduct))/((double)faDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(dialFaWoFPColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (includeChains && numberSns>1 && faDetectablesClass>0) sumCell.setCellValue(faDetectedDialWOFPsAdduct+"/"+faDetectablesAdduct);
          else sumCell.setCellValue("NA");
          sumCell = sumRow.createCell(dialFaWoFPPercentColumn);
          sumCell.setCellType(Cell.CELL_TYPE_STRING);
          sumCell.setCellStyle(centerStyle);
          if (includeChains && numberSns>1 && faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedDialWOFPsAdduct))/((double)faDetectablesAdduct), 0)+"%");
          else sumCell.setCellValue("NA");

          
          ms1Detectables += ms1DetectablesAdduct;
          faDetectables += faDetectablesAdduct;
          positionDetectables += positionDetectablesAdduct;
          ms1DetectedLDA += ms1DetectedLDAAdduct;
          faDetectedLDAWOFPs += faDetectedLDAWOFPsAdduct;
          faDetectedLDA += faDetectedLDAAdduct;
          ms1DetectedDial += ms1DetectedDialAdduct;
          faDetectedDialWOFPs += faDetectedDialWOFPsAdduct;
          faDetectedDial += faDetectedDialAdduct;
        }
        rowCountDetails++;
        row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;
        cell = row.createCell(ldaHitsColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA");
        cell = row.createCell(lbHitsColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("MS-DIAL");
                
        row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;        
        cell = row.createCell(typeColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Identified Species");
        cell = row.createCell(ldaHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        cell.setCellValue(ms1DetectedLDAClass+"/"+ms1DetectablesClass);
        cell = row.createCell(ldaPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (ms1DetectablesClass==0)
          cell.setCellValue("0%");
        else
          cell.setCellValue(Calculator.FormatNumberToString(((double)(100*ms1DetectedLDAClass))/((double)ms1DetectablesClass), 0)+"%");
        cell = row.createCell(lbHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        cell.setCellValue(ms1DetectedDialClass+"/"+ms1DetectablesClass);
        cell = row.createCell(lbPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (ms1DetectablesClass==0)
          cell.setCellValue("0%");
        else
          cell.setCellValue(Calculator.FormatNumberToString(((double)(100*ms1DetectedDialClass))/((double)ms1DetectablesClass), 0)+"%");
        
        row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;        
        cell = row.createCell(typeColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Identified chain compositions");
        cell = row.createCell(ldaHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && faDetectablesClass>0) cell.setCellValue(faDetectedLDAClass+"/"+faDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(ldaPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && faDetectablesClass>0) cell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLDAClass))/((double)faDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");
        cell = row.createCell(lbHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && faDetectablesClass>0) cell.setCellValue(faDetectedDialClass+"/"+faDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(lbPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && faDetectablesClass>0) cell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedDialClass))/((double)faDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");

        row = detailsSheet.createRow(rowCountDetails);
        rowCountDetails++;        
        cell = row.createCell(typeColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Identified chain compositions w.o. FPs");
        cell = row.createCell(ldaHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && faDetectablesClass>0) cell.setCellValue(faDetectedLDAWOFPsClass+"/"+faDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(ldaPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && faDetectablesClass>0) cell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedLDAWOFPsClass))/((double)faDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");
        cell = row.createCell(lbHitsColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && faDetectablesClass>0) cell.setCellValue(faDetectedDialWOFPsClass+"/"+faDetectablesClass);
        else cell.setCellValue("NA");
        cell = row.createCell(lbPercentColumn);
        cell.setCellType(Cell.CELL_TYPE_STRING);
        if (numberSns>1 && ms1DetectablesClass>0 && faDetectablesClass>0) cell.setCellValue(Calculator.FormatNumberToString(((double)(100*faDetectedDialWOFPsClass))/((double)faDetectablesClass), 0)+"%");
        else cell.setCellValue("NA");

        rowCountDetails++;
        rowCountDetails++;

        if (ms1DetectablesClass==0) continue;
        speciesEvidence.put(lipidClass, speciesEvidenceClass);
                
      }
      rowCountSummary++;
      
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("Total");
      
      sumCell = sumRow.createCell(ldaMs1Column);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(ms1DetectedLDA+"/"+ms1Detectables);
      sumCell = sumRow.createCell(ldaMs1PercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLDA)/((double)ms1Detectables), 0)+"%");
      sumCell = sumRow.createCell(dialMs1Column);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(ms1DetectedDial+"/"+ms1Detectables);
      sumCell = sumRow.createCell(dialMs1PercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedDial)/((double)ms1Detectables), 0)+"%");
      sumCell = sumRow.createCell(ldaFaColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLDA+"/"+faDetectables);
      sumCell = sumRow.createCell(ldaFaPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDA)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(dialFaColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedDial+"/"+faDetectables);
      sumCell = sumRow.createCell(dialFaPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedDial)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(ldaFaWoFPColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLDAWOFPs+"/"+faDetectables);
      sumCell = sumRow.createCell(ldaFaWoFPPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDAWOFPs)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(dialFaWoFPColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedDialWOFPs+"/"+faDetectables);
      sumCell = sumRow.createCell(dialFaWoFPPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedDialWOFPs)/((double)faDetectables), 0)+"%");

      rowCountSummary++;
      rowCountSummary++;
      rowCountSummary++;
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(ldaMs1Column);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA identified");
      sumCell = sumRow.createCell(dialMs1Column);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("MS-DIAL identified");
      sumCell = sumRow.createCell(ldaFaColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA FAs");
      sumCell = sumRow.createCell(dialFaColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("MS-DIAL FAs");
      sumCell = sumRow.createCell(ldaFaWoFPColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("LDA FAs no FPs");
      sumCell = sumRow.createCell(dialFaWoFPColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("MS-DIAL FAs no FPs");
      
      ms1Detectables = 0;   
      faDetectables = 0;
      positionDetectables = 0;  
      ms1DetectedLDA = 0;
      faDetectedLDAWOFPs = 0;
      faDetectedLDA = 0;
      ms1DetectedDial = 0;
      faDetectedDialWOFPs = 0;
      faDetectedDial = 0;
      for (String lipidClass : speciesEvidence.keySet()){
        sumRow = summarySheet.createRow(rowCountSummary);
        rowCountSummary++;
        sumCell = sumRow.createCell(sumTypeColumn);
        sumCell.setCellStyle(leftHeaderStyle);
        sumCell.setCellValue(lipidClass);
        
        LinkedHashMap<String,LdaDialStandardsEvidence> evidenceClass =  speciesEvidence.get(lipidClass);
        int ms1DetectablesClass = 0;
        int faDetectablesClass = 0;
        int positionDetectablesClass = 0;
        int ms1DetectedLDAClass = 0;
        int faDetectedLDAWOFPsClass = 0;
        int faDetectedLDAClass = 0;
        int posDetectedLDAClass = 0;
        int ms1DetectedLBClass = 0;
        int faDetectedLBWOFPsClass = 0;
        int faDetectedLBClass = 0;
        int posDetectedLBClass = 0;
        for (LdaDialStandardsEvidence ev: evidenceClass.values()){
          ms1DetectablesClass++;
          if (ev.isMs1DetectedLDA()) ms1DetectedLDAClass++;
          if (ev.isMs1DetectedLB()) ms1DetectedLBClass++;
          if (ev.isFaDetectable()){
            faDetectablesClass++;
            if (ev.isFaDetectedLDA()) faDetectedLDAClass++;
            if (ev.isFaDetectedLB()) faDetectedLBClass++;
            if (ev.isFaDetectedLDAWoFPs()) faDetectedLDAWOFPsClass++;
            if (ev.isFaDetectedLBWoFPs()) faDetectedLBWOFPsClass++;
          }
        }
        
        sumCell = sumRow.createCell(ldaMs1Column);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(ms1DetectedLDAClass+"/"+ms1DetectablesClass);
        sumCell = sumRow.createCell(ldaMs1PercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLDAClass)/((double)ms1DetectablesClass), 0)+"%");
        sumCell = sumRow.createCell(dialMs1Column);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(ms1DetectedLBClass+"/"+ms1DetectablesClass);
        sumCell = sumRow.createCell(dialMs1PercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLBClass)/((double)ms1DetectablesClass), 0)+"%");
        sumCell = sumRow.createCell(ldaFaColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(faDetectedLDAClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaFaPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDAClass)/((double)faDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(dialFaColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(faDetectedLBClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(dialFaPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLBClass)/((double)faDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaFaWoFPColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(faDetectedLDAWOFPsClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(ldaFaWoFPPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDAWOFPsClass)/((double)faDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(dialFaWoFPColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(faDetectedLBWOFPsClass+"/"+faDetectablesClass);
        else sumCell.setCellValue("NA");
        sumCell = sumRow.createCell(dialFaWoFPPercentColumn);
        sumCell.setCellType(Cell.CELL_TYPE_STRING);
        sumCell.setCellStyle(centerStyle);
        if (faDetectablesClass>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLBWOFPsClass)/((double)faDetectablesClass), 0)+"%");
        else sumCell.setCellValue("NA");

        ms1Detectables += ms1DetectablesClass;
        faDetectables += faDetectablesClass;
        positionDetectables += positionDetectablesClass;
        ms1DetectedLDA += ms1DetectedLDAClass;
        faDetectedLDAWOFPs += faDetectedLDAWOFPsClass;
        faDetectedLDA += faDetectedLDAClass;
        ms1DetectedDial += ms1DetectedLBClass;
        faDetectedDialWOFPs += faDetectedLBWOFPsClass;
        faDetectedDial += faDetectedLBClass;
      }
      rowCountSummary++;
      
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("Total");
      
      sumCell = sumRow.createCell(ldaMs1Column);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(ms1DetectedLDA+"/"+ms1Detectables);
      sumCell = sumRow.createCell(ldaMs1PercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedLDA)/((double)ms1Detectables), 0)+"%");
      sumCell = sumRow.createCell(dialMs1Column);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(ms1DetectedDial+"/"+ms1Detectables);
      sumCell = sumRow.createCell(dialMs1PercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*ms1DetectedDial)/((double)ms1Detectables), 0)+"%");
      sumCell = sumRow.createCell(ldaFaColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLDA+"/"+faDetectables);
      sumCell = sumRow.createCell(ldaFaPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDA)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(dialFaColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedDial+"/"+faDetectables);
      sumCell = sumRow.createCell(dialFaPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedDial)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(ldaFaWoFPColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedLDAWOFPs+"/"+faDetectables);
      sumCell = sumRow.createCell(ldaFaWoFPPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedLDAWOFPs)/((double)faDetectables), 0)+"%");
      sumCell = sumRow.createCell(dialFaWoFPColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(faDetectedDialWOFPs+"/"+faDetectables);
      sumCell = sumRow.createCell(dialFaWoFPPercentColumn);
      sumCell.setCellType(Cell.CELL_TYPE_STRING);
      sumCell.setCellStyle(centerStyle);
      if (faDetectables>0) sumCell.setCellValue(Calculator.FormatNumberToString((double)(100*faDetectedDialWOFPs)/((double)faDetectables), 0)+"%");
      
      //legend
      rowCountSummary++;
      rowCountSummary++;
      rowCountSummary++;
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(leftHeaderStyle);
      sumCell.setCellValue("Legend to the \"Details\" tab:");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(fullyCorrectStyle);
      sumCell.setCellValue("Green corresponds to a completely correct lipid species identifications");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(faCorrectStyle);
      sumCell.setCellValue("Cyan corresponds to correctly identified fatty acid (FA) position, without any false positives (FPs), and without  position assignment");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(faFoundStyle);
      sumCell.setCellValue("Dark yellow corresponds to correctly identified FA positions, but with additional FPs");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(ms1FoundStyle);
      sumCell.setCellValue("Orange corresponds to a detection of the lipid class, but without identifying the correct FAs");
      sumRow = summarySheet.createRow(rowCountSummary);
      rowCountSummary++;
      sumCell = sumRow.createCell(sumTypeColumn);
      sumCell.setCellStyle(notFoundStyle);
      sumCell.setCellValue("Red corresponds to an undetected analyte (this color is not used if no MS/MS are present)");

      summarySheet.setColumnWidth(0, 5000);
      
      resultWorkbook.write(out);
      resultWorkbook.close();
      out.close();
    }catch(Exception ex){
      ex.printStackTrace();
    } 
  }
  
  private void parseMSFinderStructure() {
    MSFinderStructureParser finderParser = new MSFinderStructureParser("C:\\Sphingolipids\\Experiment1\\Obitrap\\negative\\MS-Dial_Mix1\\export_MSFinder2\\Structure result-2101.txt");
    try {
      finderParser.parse();
//      Vector<MSDialEntry> entries = dialParser.getResults();
//      for (MSDialEntry entry : entries) {
//        System.out.println(entry.getName());
//      }
      
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  private void generateIsoLabeledMassList() {
    try{
      ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
      parser.parse();
      FALibParser faParser = new FALibParser("C:\\Development\\LipidDataAnalyzer\\fattyAcids\\fattyAcidChains.xlsx");
      faParser.parseFile();
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,FattyAcidVO>>> fas = faParser.getFattyAcids().get("n");
      Vector<FattyAcidVO> chains = new Vector<FattyAcidVO>();
      for (Integer cAtoms : fas.keySet()) {
        for (Integer dbs : fas.get(cAtoms).keySet()) {
          for (String pref : fas.get(cAtoms).get(dbs).keySet()) {
            chains.add(fas.get(cAtoms).get(dbs).get(pref));
          }
        }
      }
      //Vector<String> combis = FragmentCalculator.getAllPossibleChainCombinations(2,chains);
      
      // for PC 
      //String quantFile = "C:\\LDA-Collaborations\\Joe\\massLists\\PC_negative.xlsx";
      // for PE 
      //String quantFile = "C:\\LDA-Collaborations\\Joe\\massLists\\PE_negative.xlsx";
      // for PS 
      //String quantFile = "C:\\LDA-Collaborations\\Joe\\massLists\\PS_negative.xlsx";
      // for PI 
      //String quantFile = "C:\\LDA-Collaborations\\Joe\\massLists\\PI_negative.xlsx";
      // for PG 
      //String quantFile = "C:\\LDA-Collaborations\\Joe\\massLists\\PG_negative.xlsx";
      // for P-PE 
      //String quantFile = "C:\\LDA-Collaborations\\Joe\\massLists\\P-PE_negative.xlsx";
      // for LPC 
      //String quantFile = "C:\\LDA-Collaborations\\Joe\\massLists\\LPC_negative.xlsx";
      // for LPE 
      //String quantFile = "C:\\LDA-Collaborations\\Joe\\massLists\\LPE_negative.xlsx";
      // for LPE 
      String quantFile = "C:\\LDA-Collaborations\\Joe\\massLists\\LPS_negative.xlsx";
      
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(quantFile));
      XSSFWorkbook resultWorkbook = new XSSFWorkbook();
      XSSFCellStyle headerStyle = getHeaderStyle(resultWorkbook);
      
      //for PC mass list
      //Sheet resultSheet = resultWorkbook.createSheet("PC");
      //for PE mass list
      //Sheet resultSheet = resultWorkbook.createSheet("PE");
      //for PS mass list
      //Sheet resultSheet = resultWorkbook.createSheet("PS");
      //for PI mass list
      //Sheet resultSheet = resultWorkbook.createSheet("PI");      
      //for PG mass list
      //Sheet resultSheet = resultWorkbook.createSheet("PG");      
      //for P-PE mass list
      //Sheet resultSheet = resultWorkbook.createSheet("P-PE");      
      //for P-PE mass list
      //Sheet resultSheet = resultWorkbook.createSheet("LPC");      
      //for LPE mass list
      //Sheet resultSheet = resultWorkbook.createSheet("LPE");      
      //for LPS mass list
      Sheet resultSheet = resultWorkbook.createSheet("LPS");      

      
      int rowCount = 4;
      Row outRow = resultSheet.createRow(rowCount);
      Cell label = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("Name");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("dbs");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("C");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("H");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("O");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("P");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("N");
      label.setCellStyle(headerStyle);
      label = outRow.createCell(8,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("D");
      label.setCellStyle(headerStyle);
      //for sphingoid backbone chain library
//      label = outRow.createCell(9,HSSFCell.CELL_TYPE_STRING);
//      label.setCellValue("mass");
//      label.setCellStyle(headerStyle);
      //for all others
      label = outRow.createCell(9,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("M");
      label.setCellStyle(headerStyle);

      
      //for positive ion mode data
//      label = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);      
//      label.setCellValue("mass(form[+H] name[H])");
//      label.setCellStyle(headerStyle);
//      label = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);      
//      label.setCellValue("mass(form[-OH] name[-OH])");
//      label.setCellStyle(headerStyle);
//      label = outRow.createCell(12,HSSFCell.CELL_TYPE_STRING);      
//      label.setCellValue("mass(form[+Na] name[Na])");
//      label.setCellStyle(headerStyle);
      //for negative ion mode data
      label = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);      
      label.setCellValue("mass(form[-H] name[-H])");
      label.setCellStyle(headerStyle);
//      label = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);      
//      label.setCellValue("mass(form[+HCOO] name[HCOO])");
//      label.setCellStyle(headerStyle);
//      label = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);      
//      label.setCellValue("mass(form[-CH3] name[-CH3])");
//      label.setCellStyle(headerStyle);

//      label = outRow.createCell(12,HSSFCell.CELL_TYPE_STRING);      
//      label.setCellValue("mass(form[+Cl] name[Cl])");
//      label.setCellStyle(headerStyle);
      
      //for all except for sphingoid backbone chain library
      label = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);
      label.setCellValue("tR (min)");
      label.setCellStyle(headerStyle);
      rowCount++;
      
      Hashtable<Integer,Integer> cDbsCombi = new Hashtable<Integer,Integer>();

      
      //for PC, PE, PS, PI, PG, P-PE
//      for (int i=20; i!=49; i++) {
//        for (int j=0; j!=13; j++) {
//          if (j>4 && i<30)
//            continue;
//          else if (j>6 && i<36)
//            continue;
//          else if (j>8 && i<38)
//            continue;
//          else if (j>10 && i<40)
//            continue;
//          cDbsCombi.put(i, j);
//        }
//      }

      //for LPC, LPE, LPS
      for (int i=12; i!=27; i++) {
        for (int j=0; j!=9; j++) {
          if (j>4 && i<16)
            continue;
          else if (j>6 && i<20)
            continue;
          cDbsCombi.put(i, j);
        }
      }

      
      //for sphingoid backbone chain library
      //for (int cAtoms=10; cAtoms<21; cAtoms+=1){
      
      //for PC, PE, PS, PI, PG, P-PE
      //for (int cAtoms=20; cAtoms<49; cAtoms+=1){
      //for LPC, LPE, LPS
      for (int cAtoms=12; cAtoms<27; cAtoms+=1){
        int dbs = 0;
        int dbsMax = cDbsCombi.get(cAtoms)+1;
        while (dbs<dbsMax){
          //for PC, PE, PS, PI, PG, P-PE
          //Vector<Vector<FattyAcidVO>> filtered = getPossCombis(2,cAtoms,dbs,chains, new Vector<FattyAcidVO>());
          //for LPC, LPE
          Vector<Vector<FattyAcidVO>> filtered = getPossCombis(1,cAtoms,dbs,chains, new Vector<FattyAcidVO>());
          int maxDeut = 0;
          for (Vector<FattyAcidVO> combs : filtered) {
            int nrDeut = 0;
            for (FattyAcidVO fa : combs) {
              if (fa.getPrefix().length()>0)
                nrDeut++;
            }
            if (nrDeut>maxDeut)
              maxDeut = nrDeut;
          }
//            for (String combi:combis) {
//              System.out.println(combi);
//            }

          for (int i=0; i!=(maxDeut+1); i++){ 
        
            //for PC chain library
//            int totalC = cAtoms+8;
//            int hAtoms = totalC*2-2*dbs;
//            int oAtoms = 8;
//            int pAtoms = 1;
//            int nAtoms = 1;

            //for PE chain library
//            int totalC = cAtoms+5;
//            int hAtoms = totalC*2-2*dbs;
//            int oAtoms = 8;
//            int pAtoms = 1;
//            int nAtoms = 1;

            //for PS chain library
//            int totalC = cAtoms+9;
//            int hAtoms = totalC*2-2*dbs-3;
//            int oAtoms = 13;
//            int pAtoms = 1;
//            int nAtoms = 0;

            //for PG chain library
//            int totalC = cAtoms+6;
//            int hAtoms = totalC*2-2*dbs-1;
//            int oAtoms = 10;
//            int pAtoms = 1;
//            int nAtoms = 0;

            //for P-PE chain library
//            int totalC = cAtoms+5;
//            int hAtoms = totalC*2-2*dbs;
//            int oAtoms = 7;
//            int pAtoms = 1;
//            int nAtoms = 1;

            //for LPC chain library
//            int totalC = cAtoms+8;
//            int hAtoms = totalC*2-2*dbs+2;
//            int oAtoms = 7;
//            int pAtoms = 1;
//            int nAtoms = 1;

//            //for LPE chain library
//            int totalC = cAtoms+5;
//            int hAtoms = totalC*2-2*dbs+2;
//            int oAtoms = 7;
//            int pAtoms = 1;
//            int nAtoms = 1;

            //for LPS chain library
            int totalC = cAtoms+6;
            int hAtoms = totalC*2-2*dbs;
            int oAtoms = 9;
            int pAtoms = 1;
            int nAtoms = 1;

            
            int dAtoms = 0;
            String name = "";
            for (int j=0; j!=i; j++) {
              name+="D";
              hAtoms -= 11;
              dAtoms += 11;
            }

            
            
            double massNeutral = parser.getElementDetails("C").getMonoMass()*((double)totalC)+parser.getElementDetails("H").getMonoMass()*((double)hAtoms)+
                parser.getElementDetails("O").getMonoMass()*((double)oAtoms)+parser.getElementDetails("P").getMonoMass()*((double)pAtoms)+
                parser.getElementDetails("N").getMonoMass()*((double)nAtoms)+parser.getElementDetails("D").getMonoMass()*((double)dAtoms);
            //double massCharged = massNeutral+parser.getElementDetails("H").getMonoMass()+parser.getElementDetails("C").getMonoMass()+2d*parser.getElementDetails("O").getMonoMass();
          
            //for positive ion mode
            //protonated
//            double massH = massNeutral+parser.getElementDetails("h").getMonoMass();
//            //water loss
//            double massOH = massNeutral+parser.getElementDetails("h").getMonoMass()-2*parser.getElementDetails("H").getMonoMass()-parser.getElementDetails("O").getMonoMass();
//            //sodiated
//            double massNa = massNeutral+parser.getElementDetails("na").getMonoMass();

          //for negative ion mode
          //deprotonated
          double massH = massNeutral-parser.getElementDetails("h").getMonoMass();
          //formate
//          double massHCOO = massNeutral-parser.getElementDetails("h").getMonoMass()+2*parser.getElementDetails("H").getMonoMass()+parser.getElementDetails("C").getMonoMass()+2*parser.getElementDetails("O").getMonoMass();  
//          //methyl loss
//          double massCH3 = massNeutral-parser.getElementDetails("h").getMonoMass()-2*parser.getElementDetails("H").getMonoMass()-parser.getElementDetails("C").getMonoMass();

            
            outRow = resultSheet.createRow(rowCount);
            label = outRow.createCell(0,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(name+cAtoms);
            label = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
            label.setCellValue(":");       
            label = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(dbs);
            label = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(totalC);
            label = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(hAtoms);
            label = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(oAtoms);
            label = outRow.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(pAtoms);
            label = outRow.createCell(7,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(nAtoms);
            label = outRow.createCell(8,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(dAtoms);
            label = outRow.createCell(9,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(massNeutral);
          
          //positive
//          label = outRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
//          label.setCellValue(massH);
//          label = outRow.createCell(11,HSSFCell.CELL_TYPE_NUMERIC);
//          label.setCellValue(massOH);
//          label = outRow.createCell(12,HSSFCell.CELL_TYPE_NUMERIC);
//          label.setCellValue(massNa);
          //negative
            label = outRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
            label.setCellValue(massH);
//            label = outRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
//            label.setCellValue(massHCOO);
//            label = outRow.createCell(11,HSSFCell.CELL_TYPE_NUMERIC);
//            label.setCellValue(massCH3);
            
//            label = outRow.createCell(12,HSSFCell.CELL_TYPE_NUMERIC);
//            label.setCellValue(massCl);
            rowCount++;
          }
          dbs++;
        }
      }
      
      resultWorkbook.write(out);
      out.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }

  }
  private Vector<Vector<FattyAcidVO>> getPossCombis(int chains, int cAtoms, int dbs, Vector<FattyAcidVO> fas, Vector<FattyAcidVO> added){
    Vector<Vector<FattyAcidVO>> combis = new Vector<Vector<FattyAcidVO>>();
    for (FattyAcidVO fa : fas) {
      Vector<FattyAcidVO> toAdd = new Vector<FattyAcidVO>(added);
      toAdd.add(fa);
      if (chains>1) {
        Vector<Vector<FattyAcidVO>> combis2 = getPossCombis(chains-1, cAtoms, dbs, fas, toAdd);
        if (combis2.size()>0) {
          combis.addAll(combis2);
        }
      }else{
        //System.out.println(toAdd.size());
        if (checkTotalCDbsValid(toAdd,cAtoms,dbs)) {
          String comb = "";
          for (FattyAcidVO fax : toAdd)
            comb += fax.getCarbonDbsId()+"_";
          //System.out.println(cAtoms+":"+dbs+": "+comb);
          combis.add(toAdd);
        }
      }
    }
    return combis;
  }
  
  private boolean checkTotalCDbsValid(Vector<FattyAcidVO> toAdd, int cAtoms, int dbs) {
    int totC = 0;
    int totD = 0;
    for (FattyAcidVO fa : toAdd) {
      totC += fa.getcAtoms();
      totD += fa.getDoubleBonds();
    }
    return (totC==cAtoms && totD==dbs);
  }
  
  
  private void getValidSphingoOrbitrapCIDSpeciesPositive(LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses,
      Hashtable<String,LipidClassInfoVO> lipidClassInfo, LinkedHashMap<String,Boolean> adducts){
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("SphBase", BrainSpecies.getSphBaseSpeciesOrbitrap());
    adducts.put("H", false);
    lipidClassInfo.put("SphBase", new LipidClassInfoVO(1,true,0.4d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("Cer", BrainSpecies.getCerSpeciesOrbitrap());
    adducts.put("H", true);
//    adducts.put("Na", false);
    lipidClassInfo.put("Cer", new LipidClassInfoVO(2,true,0.4d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("Cer1P", BrainSpecies.getCer1PSpeciesOrbitrap());
    adducts.put("H", true);
//    adducts.put("Na", false);
    lipidClassInfo.put("Cer1P", new LipidClassInfoVO(2,false,0.4d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("HexCer", BrainSpecies.getHexCerSpeciesOrbitrap());
    adducts.put("H", true);
//    adducts.put("Na", false);
    lipidClassInfo.put("HexCer", new LipidClassInfoVO(2,true,0.4d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("SM", BrainSpecies.getSMSpeciesOrbitrap());
    adducts.put("H", false);
    adducts.put("Na", false);
    lipidClassInfo.put("SM", new LipidClassInfoVO(1,true,0.4d,adducts));    
  }
  
  private void getValidSphingoOrbitrapCIDSpeciesNegative(LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses,
      Hashtable<String,LipidClassInfoVO> lipidClassInfo, LinkedHashMap<String,Boolean> adducts){
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("Cer", BrainSpecies.getCerSpeciesOrbitrap());
    adducts.put("HCOO", false);
    adducts.put("-H", true);
    lipidClassInfo.put("Cer", new LipidClassInfoVO(2,true,0.4d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("Cer1P", BrainSpecies.getCer1PSpeciesOrbitrap());
    adducts.put("-H", true);
    lipidClassInfo.put("Cer1P", new LipidClassInfoVO(2,false,0.4d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("HexCer", BrainSpecies.getHexCerSpeciesOrbitrap());
    adducts.put("HCOO", false);
    adducts.put("-H", true);
    lipidClassInfo.put("HexCer", new LipidClassInfoVO(2,true,0.4d,adducts));
    adducts = new LinkedHashMap<String,Boolean>();
    lipidClasses.put("SM", BrainSpecies.getSMSpeciesOrbitrap());
    adducts.put("HCOO", false);
    lipidClassInfo.put("SM", new LipidClassInfoVO(2,true,0.4d,adducts));    
  }

  
  private void compareLDAMSDialNaturalProbesPositive(){
    //the first key is the lipid class, the second key the ms1 species name, the third key the structural identification
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>>();
    Hashtable<String,LipidClassInfoVO> lipidClassInfo = new Hashtable<String,LipidClassInfoVO>();
    LinkedHashMap<String,Boolean> adducts = new LinkedHashMap<String,Boolean>();
    //TODO: the implementation of this method is not complete - only a test case
    this.getValidSphingoOrbitrapCIDSpeciesPositive(lipidClasses,lipidClassInfo,adducts);
    for (LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> correctAnalytes: lipidClasses.values()) {
      for (LinkedHashMap<String,ReferenceInfoVO> classCorrect : correctAnalytes.values()) {
        for (ReferenceInfoVO info : classCorrect.values()) {
          for (int i=0; i!=info.getCorrectRts().length; i++) {
            info.getCorrectRts()[i] = info.getCorrectRts()[i]-0.05d;        }
        }
      }
    }
    ////this.getValid4000QTRAPSpeciesNegative(lipidClasses,lipidClassInfo,adducts);

    String chromFile = "C:\\Sphingolipids\\Brain\\Orbitrap\\positive\\Brain_pos_5.chrom";
    String quantFile = "C:\\Sphingolipids\\Brain\\massLists\\Orbitrap\\sphingos_positive.xlsx";
    String ldaFile = "C:\\Sphingolipids\\Brain\\Orbitrap\\positive\\Brain_pos_5_sphingos_positive.xlsx";
    String msDialFile = "C:\\Sphingolipids\\Brain\\Orbitrap\\MS-Dial_positive\\export\\Brain_pos_5.txt";
    String outputFile = "C:\\Sphingolipids\\Brain\\Orbitrap\\MS-Dial_positive\\Brain_pos_5_comp_generated.xlsx";

    performMSDialComparisonOfNaturalProbes(lipidClasses,lipidClassInfo,chromFile,quantFile,ldaFile,msDialFile,outputFile);    
  }
  
  private void compareLDAMSDialNaturalProbesNegative(){
    //the first key is the lipid class, the second key the ms1 species name, the third key the structural identification
    LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses = new LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>>();
    Hashtable<String,LipidClassInfoVO> lipidClassInfo = new Hashtable<String,LipidClassInfoVO>();
    LinkedHashMap<String,Boolean> adducts = new LinkedHashMap<String,Boolean>();
    
    //TODO: the implementation of this method is not complete - only a test case
    this.getValidSphingoOrbitrapCIDSpeciesNegative(lipidClasses,lipidClassInfo,adducts);
    ////this.getValid4000QTRAPSpeciesNegative(lipidClasses,lipidClassInfo,adducts);

    String chromFile = "C:\\Sphingolipids\\Brain\\Orbitrap\\negative\\Brain_neg_1.chrom";
    String quantFile = "C:\\Sphingolipids\\Brain\\massLists\\Orbitrap\\sphingos_negative.xlsx";
    String ldaFile = "C:\\Sphingolipids\\Brain\\Orbitrap\\negative\\Brain_neg_1_sphingos_negative.xlsx";
    String msDialFile = "C:\\Sphingolipids\\Brain\\Orbitrap\\MS-Dial_negative\\export\\Brain_neg_1.txt";
    String outputFile = "C:\\\\Sphingolipids\\\\Brain\\\\Orbitrap\\\\MS-Dial_negative\\Brain_neg_1_comp_generated.xlsx";
//    String chromFile = "D:\\BiologicalExperiment\\QTRAP\\negative\\Data20151002_QTrap_Liver-025_QTrap_Liver1-1_neg.chrom";
//    String quantFile = "D:\\BiologicalExperiment\\massLists\\negative\\QTRAP\\negative.xlsx";
//    String ldaFile = "D:\\BiologicalExperiment\\QTRAP\\negative\\Data20151002_QTrap_Liver-025_QTrap_Liver1-1_neg_negative.xlsx";
//    String lbFile = "D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\negative\\output\\Data20151002_QTrap_Liver-025_QTrap_Liver1-1_neg_MF450.mgf.tsv";
//    String outputFile = "D:\\BiologicalExperiment\\LipidBlast\\4000QTRAP\\negative\\Data20151002_QTrap_Liver-025_QTrap_Liver1-1_neg_MF450_comp_generated.xlsx";

    performMSDialComparisonOfNaturalProbes(lipidClasses,lipidClassInfo,chromFile,quantFile,ldaFile,msDialFile,outputFile);
  }

  private void  performMSDialComparisonOfNaturalProbes(LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>>> lipidClasses, 
      Hashtable<String,LipidClassInfoVO> lipidClassInfo, String chromFile, String quantFile, String ldaFile, String msDialFile, String outputFile){
    try{
      Vector quantValues = QuantificationThread.parseQuantExcelFile(quantFile,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f, true);
      ////Vector quantValues = null;
      String[] chromPaths = StringUtils.getChromFilePaths(chromFile);   
      LipidomicsAnalyzer analyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false);

      Hashtable<String,Vector<String>> analyteSequence = (Hashtable<String,Vector<String>>)quantValues.get(1);
      Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>> quantVOs = (Hashtable<String,Hashtable<String,Hashtable<String,QuantVO>>>) quantValues.get(3);
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outputFile));
      Workbook resultWorkbook = new XSSFWorkbook();
      CellStyle fullyCorrectStyle = getFullyCorrectStyle(resultWorkbook);
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      CellStyle leftHeaderStyle = getLeftHeaderStyle(resultWorkbook);
      CellStyle boldStyle = getBoldStyle(resultWorkbook);
      CellStyle notFoundStyle = getNotFoundStyle(resultWorkbook);
      CellStyle ambigStyle = getFACorrectStyle(resultWorkbook);
      CellStyle onlyMS1Style = getMS1FoundStyle(resultWorkbook);
      CellStyle percentageStyle = getPercentageStyle(resultWorkbook);
      leftHeaderStyle.cloneStyleFrom(headerStyle);
      leftHeaderStyle.setAlignment(CellStyle.ALIGN_LEFT);
      Hashtable<String,Vector<LipidParameterSet>> resultsLDA = LDAResultReader.readResultFile(ldaFile, new Hashtable<String,Boolean>()).getIdentifications();
      MSDialTxtParser msdialParser = new MSDialTxtParser(msDialFile);
      msdialParser.parse();
      //Hashtable<String,Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>> resultsLB = lBlastParser.getResults_();
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> resultsDial = msdialParser.getStructuredResults();
      Hashtable<String,Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>> dialMS1Only = msdialParser.getStructuredMS1Only();
      
      //check which DIAL classes were not found
      for (String dialClass : resultsDial.keySet()) {
        if (!lipidClasses.containsKey(dialClass))
          System.out.println("!!!!! The following MS-DIAL class is not supported: "+dialClass);
      }
      
      //remove DIAL hits with four or more hydroxygens
      for (String lClass : resultsDial.keySet()) {
        Vector<String> toRemove = new Vector<String>();
        int[] cAndDbs;
        for (String analyte : resultsDial.get(lClass).keySet()) {
          if ((lClass.equalsIgnoreCase("Cer") || lClass.equalsIgnoreCase("Cer1P") || lClass.equalsIgnoreCase("HexCer")) && analyte.startsWith("q"))
            toRemove.add(analyte);
          else if (lClass.equalsIgnoreCase("Cer") || lClass.equalsIgnoreCase("Cer1P") || lClass.equalsIgnoreCase("HexCer") || lClass.equalsIgnoreCase("SM")){
            cAndDbs = StaticUtils.parseCAndDbsFromChainId(analyte.substring(1));
            if (cAndDbs[0]<26 && cAndDbs[1]>1)
              toRemove.add(analyte);
            else if (cAndDbs[0]<28 && cAndDbs[1]>2)
              toRemove.add(analyte);
            else if (cAndDbs[0]<32 && cAndDbs[1]>3)
              toRemove.add(analyte);
            else if (cAndDbs[0]<34 && cAndDbs[1]>4)
              toRemove.add(analyte);
            else if (cAndDbs[0]<36 && cAndDbs[1]>5)
              toRemove.add(analyte);
            else if (cAndDbs[0]<20 || cAndDbs[0]>48 || cAndDbs[1]>7)
              toRemove.add(analyte);
          }
        }
        for (String analyte : toRemove) {
          resultsDial.get(lClass).remove(analyte);
        }
      }
      
      for (String lClass : lipidClasses.keySet()){
        LinkedHashMap<String,LinkedHashMap<String,ReferenceInfoVO>> species = lipidClasses.get(lClass);
        LipidClassInfoVO classInfo = lipidClassInfo.get(lClass);
        Sheet sheet = resultWorkbook.createSheet(lClass);
        int rowCount = 0;
        Row row =  sheet.createRow(rowCount);
        int classColumn = 0;
        int width = (14+1)*256;
        int ms1Column = 1;
        sheet.setColumnWidth(ms1Column, width);
        int adductColumn = 2;
        width = (12+1)*256;
        int ldaIdentification = 3;
        sheet.setColumnWidth(ldaIdentification, width);
        int lbIdentification = 4;
        width = (12+1)*256;
        sheet.setColumnWidth(lbIdentification, width);
        int ldaMS1CodeColumn = 5;
        width = (12+1)*256;
        sheet.setColumnWidth(ldaMS1CodeColumn, width);
        int lbMS1CodeColumn = 6;
        width = (12+1)*256;
        sheet.setColumnWidth(lbMS1CodeColumn, width);
        int structureColumn = 7;
        int ldaIdentificationMS2 = 8;
        int lbIdentificationMS2 = 9;
        int ldaMS2CodeColumn = 10;
        int lbMS2CodeColumn = 11;
        int amountMS2Columns = 5;
        int rtColumn = 7;
        int rtLDAColumn = 8;
        int rtLBColumn = 9;
        int probLBColumn = 10;
        int mzColumn = 11;
        int commentColumn = 12;
        int ldaMS1IgnoreColumn = 13;
        int lbMS1IgnoreColumn = 14;
        int ldaMS2IgnoreColumn = 15;
        int lbMS2IgnoreColumn = 16;
        
        int ldaIdentCode;
        int dialIdentCode;
       
        int ldaMS2IdentCode;
        int dialMS2IdentCode;
        
        rowCount++;
        Cell cell = row.createCell(classColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Class");
        cell = row.createCell(ms1Column);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Name");
        cell = row.createCell(adductColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Adduct");
        cell = row.createCell(ldaIdentification);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA");
        cell = row.createCell(lbIdentification);
        width = (16+1)*256;
        sheet.setColumnWidth(lbIdentification, width);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("MS-Dial");
        cell = row.createCell(ldaMS1CodeColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDA-Code");
        cell = row.createCell(lbMS1CodeColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("MS-Dial-Code");
        if (classInfo.getSns()>1){
          width = (20+1)*256;
          sheet.setColumnWidth(structureColumn, width);
          cell = row.createCell(structureColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("Structure");
          width = (12+1)*256;
          sheet.setColumnWidth(ldaIdentificationMS2, width);
          cell = row.createCell(ldaIdentificationMS2);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("LDA");
          width = (22+1)*256;
          sheet.setColumnWidth(lbIdentificationMS2, width);
          cell = row.createCell(lbIdentificationMS2);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("MS-Dial");
          width = (12+1)*256;
          sheet.setColumnWidth(ldaMS2CodeColumn, width);
          cell = row.createCell(ldaMS2CodeColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("LDA-Code");
          width = (12+1)*256;
          sheet.setColumnWidth(lbMS2CodeColumn, width);
          cell = row.createCell(lbMS2CodeColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("MS-Dial-Code");
         
          rtColumn += amountMS2Columns;
          rtLDAColumn += amountMS2Columns;
          rtLBColumn += amountMS2Columns;
          probLBColumn += amountMS2Columns;
          mzColumn += amountMS2Columns;
          commentColumn += amountMS2Columns;
          ldaMS1IgnoreColumn += amountMS2Columns;
          lbMS1IgnoreColumn += amountMS2Columns;
          ldaMS2IgnoreColumn += amountMS2Columns;
          lbMS2IgnoreColumn += amountMS2Columns;
          
          width = (16+1)*256;
          sheet.setColumnWidth(ldaMS2IgnoreColumn, width);
          cell = row.createCell(ldaMS2IgnoreColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("LDAMS2-Ignore");
          width = (16+1)*256;
          sheet.setColumnWidth(lbMS2IgnoreColumn, width);
          cell = row.createCell(lbMS2IgnoreColumn);
          cell.setCellStyle(headerStyle);
          cell.setCellValue("MSDMS2-Ignore");
        }
        width = (9+1)*256;
        sheet.setColumnWidth(rtColumn, width);
        cell = row.createCell(rtColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("RT-ideal");
        width = (9+1)*256;
        sheet.setColumnWidth(rtLDAColumn, width);
        cell = row.createCell(rtLDAColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("RT-LDA");
        width = (9+1)*256;
        sheet.setColumnWidth(rtLBColumn, width);
        cell = row.createCell(rtLBColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("RT-MS-Dial");
        width = (9+1)*256;
        sheet.setColumnWidth(probLBColumn, width);
        cell = row.createCell(probLBColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("MSD-Prob");
        width = (26+1)*256;
        sheet.setColumnWidth(commentColumn, width);
        cell = row.createCell(mzColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("m/z");
        width = (30+1)*256;
        sheet.setColumnWidth(commentColumn, width);
        cell = row.createCell(commentColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("Comment");      
        width = (16+1)*256;
        sheet.setColumnWidth(ldaMS1IgnoreColumn, width);
        cell = row.createCell(ldaMS1IgnoreColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("LDAMS1-Ignore");
        width = (16+1)*256;
        sheet.setColumnWidth(lbMS1IgnoreColumn, width);
        cell = row.createCell(lbMS1IgnoreColumn);
        cell.setCellStyle(headerStyle);
        cell.setCellValue("MSDMS1-Ignore");


        int finalRowNumberAnalytes = 0;
        Hashtable<String,Hashtable<String,QuantVO>> quantsOfClass = quantVOs.get(lClass);
        Vector<LipidParameterSet> ldaAnalytes = resultsLDA.get(lClass);
        Hashtable<String,Hashtable<String,Vector<MSDialEntry>>> dialAnalytes = new Hashtable<String,Hashtable<String,Vector<MSDialEntry>>>();
        if (resultsDial.containsKey(lClass)) dialAnalytes = resultsDial.get(lClass);
        Vector<String> dialsNotInQuantFile = new Vector<String>(); 
        //check which DIAL analytes were not found
        for (String analyte : dialAnalytes.keySet()) {
          boolean found = false;
          for (String inQuant : analyteSequence.get(lClass)) {
            if (analyte.equalsIgnoreCase(inQuant)) {
              found = true;
              break;
            }
          }
          if (!found)
            dialsNotInQuantFile.add(analyte);
        }
        
        Vector<String> allAnalytes = new Vector<String>(analyteSequence.get(lClass));
        //put the missing MS-Dial hits at the correct position in the allAnalytes table
        for (String analyte : dialsNotInQuantFile) {
          //int[0] nr of C atoms; int[1] nr of dbs; int[2] nr of OH
          int[] cAtomsDbsAndOh = decodeAnalyte(analyte);
          int addPosition = 0;
          for (String other : allAnalytes) {
            int[] otherCAtomsDbsAndOh = decodeAnalyte(other);
            if (cAtomsDbsAndOh[0]<otherCAtomsDbsAndOh[0])
              break;
            else if (cAtomsDbsAndOh[0]==otherCAtomsDbsAndOh[0]) {
              if (cAtomsDbsAndOh[1]<otherCAtomsDbsAndOh[1])
                break;
              else if (cAtomsDbsAndOh[1]==otherCAtomsDbsAndOh[1]&&cAtomsDbsAndOh[2]<otherCAtomsDbsAndOh[2])
                break;
            }
            addPosition++;
          }
          allAnalytes.add(addPosition,analyte);
        }
        for (String analyte : allAnalytes){
          if (analyte.startsWith("IS")) continue;
          Hashtable<String,QuantVO> quantsOfAnalyte = quantsOfClass.get(analyte);
          Hashtable<String,Vector<MSDialEntry>> dialsOfMod = new Hashtable<String,Vector<MSDialEntry>>();
          if (dialAnalytes.containsKey(analyte)) dialsOfMod = dialAnalytes.get(analyte);
          LinkedHashMap<String,ReferenceInfoVO> molecularSpecies = null;
          if (species.containsKey(analyte)) molecularSpecies = species.get(analyte);
          boolean doNotUse = true;
          if (molecularSpecies==null || molecularSpecies.size()==0)
            doNotUse = false;
          else {
            for (ReferenceInfoVO info : molecularSpecies.values()) {
              if (info.useInEvaluation()) doNotUse = false;
            }
          }
          Vector<LdaLBLASTCompareVO> comparesOfOneSpecies = new Vector<LdaLBLASTCompareVO>();
          
         for (String adduct : classInfo.getAdducts().keySet()){
            //if (!quantsOfAnalyte.containsKey(adduct)) continue;
           Vector<LipidParameterSet> ldaAnalyte = getLDAAnalytes(ldaAnalytes, analyte, adduct, true);
           Vector<MSDialEntry> ident = dialsOfMod.get(adduct);
           if (ident!=null || ldaAnalyte.size()>0){
             float mz = 0f;
             if (quantsOfAnalyte!=null && quantsOfAnalyte.containsKey(adduct))
               mz = (float)quantsOfAnalyte.get(adduct).getAnalyteMass();
             else if (ldaAnalyte.size()>0)
               mz = ldaAnalyte.get(0).Mz[0];
             else if (ident!=null && ident.size()>0)
               mz = (float)ident.get(0).getMz();
             Vector<MSDialEntry> dialMS1 = null;
             if (dialMS1Only.containsKey(lClass) && dialMS1Only.get(lClass).containsKey(analyte) && dialMS1Only.get(lClass).get(analyte).containsKey(adduct))
               dialMS1 = dialMS1Only.get(lClass).get(analyte).get(adduct);
             LdaLBLASTCompareVO ms1Compare = getLdaMSDialSpeciesComparision(analyte,classInfo,molecularSpecies,ldaAnalyte,ident,dialMS1,analyzer,mz,lClass,classInfo.getAdducts().get(adduct));
             if (ms1Compare==null) continue;
             ms1Compare.setAdduct(adduct);
             comparesOfOneSpecies.add(ms1Compare);
           }
          }
          setIgnoresForSpeciesLevelEvaluation(comparesOfOneSpecies,classInfo.getSns()>1);
          for (LdaLBLASTCompareVO ms1Compare : comparesOfOneSpecies){
            row =  sheet.createRow(rowCount);
            rowCount++;
            cell = row.createCell(classColumn);
            cell.setCellType(Cell.CELL_TYPE_STRING);
            cell.setCellValue(lClass);
            cell = row.createCell(ms1Column);
            cell.setCellType(Cell.CELL_TYPE_STRING);
            if (ms1Compare.isFp()) cell.setCellStyle(notFoundStyle);
            else cell.setCellStyle(fullyCorrectStyle );
            cell.setCellValue(ms1Compare.getCorrectName()+(doNotUse?"_noDial":""));
            cell = row.createCell(adductColumn);
            cell.setCellType(Cell.CELL_TYPE_STRING);
            cell.setCellValue(ms1Compare.getAdduct());
              
            CellStyle ldaMS1Style = getStyleBasedOnEvidence(ms1Compare.getLdaIdentCode(),ms1Compare.isFp(),fullyCorrectStyle,
                ambigStyle,onlyMS1Style,notFoundStyle);
            CellStyle lbMS1Style = getStyleBasedOnEvidence(ms1Compare.getLbIdentCode(),ms1Compare.isFp(),fullyCorrectStyle,
                ambigStyle,onlyMS1Style,notFoundStyle);
              
              
              cell = row.createCell(ldaIdentification);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              if (ldaMS1Style!=null) cell.setCellStyle(ldaMS1Style);
              cell.setCellValue(ms1Compare.getLdaName());
              cell = row.createCell(lbIdentification);
              cell.setCellType(Cell.CELL_TYPE_STRING);
              if (lbMS1Style!=null) cell.setCellStyle(lbMS1Style);
              cell.setCellValue(ms1Compare.getLbName());
              cell = row.createCell(ldaMS1CodeColumn);
              cell.setCellType(Cell.CELL_TYPE_NUMERIC);
              if (ldaMS1Style!=null) cell.setCellStyle(ldaMS1Style);
              ldaIdentCode = ms1Compare.getLdaIdentCode();
              if (doNotUse)
                ldaIdentCode = 0;
              cell.setCellValue(ldaIdentCode);
              cell = row.createCell(lbMS1CodeColumn);
              cell.setCellType(Cell.CELL_TYPE_NUMERIC);
              if (lbMS1Style!=null) cell.setCellStyle(lbMS1Style);
              cell.setCellValue(ms1Compare.getLbIdentCode());
              if (ms1Compare.isIgnoreLDA()){
                cell = row.createCell(ldaMS1IgnoreColumn);
                cell.setCellValue(true);
              }
              if (ms1Compare.isIgnoreLB()){
                cell = row.createCell(lbMS1IgnoreColumn);
                cell.setCellValue(true);
              }
              

//              if (analyte.equalsIgnoreCase("t40:0")) {
//                System.out.println("t40:0: "+ms1Compare.getRt());
//              }

              if (classInfo.getSns()<2 || ms1Compare.getMs2Evidence().size()==0 || !classInfo.getAdducts().get(ms1Compare.getAdduct()) || doNotUse){
                cell = row.createCell(rtColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(ms1Compare.getRt());
                cell = row.createCell(rtLDAColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(ms1Compare.getLdaRts());
                cell = row.createCell(rtLBColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(ms1Compare.getLbRts());
                if (ms1Compare.getLbIdentCode()>0 || ms1Compare.getLbIdentCode()<-1){
                  cell = row.createCell(probLBColumn);
                  cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(ms1Compare.getLbProb());
                }
                cell = row.createCell(mzColumn);
                cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                cell.setCellValue(ms1Compare.getMz());                
                if (ms1Compare.isLdaMs1Only() || ms1Compare.areFasInOtherCombination() || doNotUse){
                  cell = row.createCell(commentColumn);
                  cell.setCellType(Cell.CELL_TYPE_STRING);
                  if (ms1Compare.isLdaMs1Only())
                    cell.setCellValue("not counted: only MS1 identification");
                  else if (ms1Compare.areFasInOtherCombination())
                    cell.setCellValue("not counted: only possible as additional combination");
                  else if (doNotUse) {
                    cell.setCellValue("not counted: MS-DIAL cannot identify monohydroxylated species");
                  }
                }

              }
              
              boolean chainInfo = false;
              if (classInfo.getSns()>1 && !classInfo.getAdducts().get(ms1Compare.getAdduct())) {
                for (LdaLBLASTCompareVO comp : ms1Compare.getMs2Evidence()) {
                  if (comp.getLdaIdentCode()>0 || comp.getLbIdentCode()>0)
                    chainInfo = true;
                }
              }
              for (int i=0; i!=ms1Compare.getMs2Evidence().size(); i++){
                LdaLBLASTCompareVO comp = ms1Compare.getMs2Evidence().get(i);
                if (i!=0){
                  row =  sheet.createRow(rowCount);
                  rowCount++;
                }
                cell = row.createCell(structureColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                if (comp.isFp()) cell.setCellStyle(notFoundStyle);
                else cell.setCellStyle(fullyCorrectStyle );
                cell.setCellValue(comp.getCorrectName());

                CellStyle ldaMS2Style = getStyleBasedOnEvidence(comp.getLdaIdentCode(),comp.isFp(),fullyCorrectStyle,
                    ambigStyle,onlyMS1Style,notFoundStyle);
                CellStyle lbMS2Style = getStyleBasedOnEvidence(comp.getLbIdentCode(),comp.isFp(),fullyCorrectStyle,
                    ambigStyle,onlyMS1Style,notFoundStyle);
                cell = row.createCell(ldaIdentificationMS2);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                if (ldaMS1Style!=null) cell.setCellStyle(ldaMS2Style);
                cell.setCellValue(comp.getLdaName());
                cell = row.createCell(lbIdentificationMS2);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                if (lbMS1Style!=null) cell.setCellStyle(lbMS2Style);
                cell.setCellValue(comp.getLbName());
                cell = row.createCell(ldaMS2CodeColumn);
                cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                if (ldaMS2Style!=null) cell.setCellStyle(ldaMS2Style);
                ldaMS2IdentCode = comp.getLdaIdentCode();
                if (!classInfo.getAdducts().get(ms1Compare.getAdduct()) && ldaMS2IdentCode!=0) {
                    ldaMS2IdentCode = 0;
                }
                cell.setCellValue(ldaMS2IdentCode);
                cell = row.createCell(lbMS2CodeColumn);
                cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                if (lbMS1Style!=null) cell.setCellStyle(lbMS2Style);
                dialMS2IdentCode = comp.getLbIdentCode();
                if (!classInfo.getAdducts().get(ms1Compare.getAdduct()) && dialMS2IdentCode!=0) {
                    dialMS2IdentCode = 0;
                }
                cell.setCellValue(dialMS2IdentCode);
                if (classInfo.getAdducts().get(ms1Compare.getAdduct()) || chainInfo) {
                  cell = row.createCell(rtColumn);
                  cell.setCellType(Cell.CELL_TYPE_STRING);
                  cell.setCellValue(comp.getRt());
                }
                if (ldaMS2IdentCode!=0) {
                  cell = row.createCell(rtLDAColumn);
                  cell.setCellType(Cell.CELL_TYPE_STRING);
                  cell.setCellValue(comp.getLdaRts());
                }
                cell = row.createCell(rtLBColumn);
                cell.setCellType(Cell.CELL_TYPE_STRING);
                cell.setCellValue(comp.getLbRts());
                if (comp.getLbIdentCode()>0 || comp.getLbIdentCode()<-1){
                  cell = row.createCell(probLBColumn);
                  cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(comp.getLbProb());
                }
                cell = row.createCell(mzColumn);
                cell.setCellType(Cell.CELL_TYPE_NUMERIC);
                cell.setCellValue(comp.getMz());
                if (comp.isIgnoreLDA()){
                  cell = row.createCell(ldaMS2IgnoreColumn);
                  cell.setCellValue(true);
                }
                if (comp.isIgnoreLB()){
                  cell = row.createCell(lbMS2IgnoreColumn);
                  cell.setCellValue(true);
                }
                if (!classInfo.getAdducts().get(ms1Compare.getAdduct())) {
                  cell = row.createCell(commentColumn);
                  cell.setCellType(Cell.CELL_TYPE_STRING);
                  cell.setCellValue("structure not counted: MS-DIAL cannot identify structures from MS3");
                }
              }
          }
        }
        finalRowNumberAnalytes = rowCount;

        rowCount++;
        rowCount++;        
        // now starts the summary section of the details tab
        int legendColumn = 0;
        int adductResultColumn = 5;
        int legend2Column = 7;
        int classResultColumn = 10;
        
        String ldaIgnore = "N";
        String lbIgnore = "O";
        if (classInfo.getSns()>1){
          ldaIgnore = "S";
          lbIgnore = "T";          
        }
        String ldaMS2Ignore = "U";
        String dialMS2Ignore = "V";          
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Species evaluation");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Species evaluation - adduct insensitive");
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        int totalSumRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified species:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(B2:B"+finalRowNumberAnalytes+",\"<>*_*\",B2:B"+finalRowNumberAnalytes+",\"<>\")");       
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified species:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(B2:B"+finalRowNumberAnalytes+",\"<>*_*\",B2:B"+finalRowNumberAnalytes+",\"<>\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        int ldaCorrectRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified LDA species:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\">0\")");       
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified LDA species:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\">0\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        int lbCorrectRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified MS-DIAL species:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\">0\")");       
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Total number of identified MS-DIAL species:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\">0\","+lbIgnore+"2:"+lbIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        int ldaFPRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False positives LDA:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"<>-1\",F2:F"+finalRowNumberAnalytes+",\"<>0\",F2:F"+finalRowNumberAnalytes+",\"<2\")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False positives LDA:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"<>-1\",F2:F"+finalRowNumberAnalytes+",\"<>0\",F2:F"+finalRowNumberAnalytes+",\"<2\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        int lbFPRow = rowCount;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False positives MS-DIAL:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"<>-1\",G2:G"+finalRowNumberAnalytes+",\"<>0\",G2:G"+finalRowNumberAnalytes+",\"<2\")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False positives MS-DIAL:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"<>-1\",G2:G"+finalRowNumberAnalytes+",\"<>0\",G2:G"+finalRowNumberAnalytes+",\"<2\","+lbIgnore+"2:"+lbIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False negatives LDA:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"=-1\")+COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"=-3\")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False negatives LDA:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"=-1\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")+COUNTIFS(F2:F"+finalRowNumberAnalytes+",\"=-3\","+ldaIgnore+"2:"+ldaIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False negatives MS-DIAL:");
        cell = row.createCell(adductResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"=-1\")+COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"=-3\")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("False negatives MS-DIAL:");
        cell = row.createCell(classResultColumn);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"=-1\","+lbIgnore+"2:"+lbIgnore+finalRowNumberAnalytes+",\"<>TRUE\")+COUNTIFS(G2:G"+finalRowNumberAnalytes+",\"=-3\","+lbIgnore+"2:"+lbIgnore+finalRowNumberAnalytes+",\"<>TRUE\")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Sensitivity LDA:");
        cell = row.createCell(adductResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("F"+ldaCorrectRow+"/F"+totalSumRow);
        cell = row.createCell(legend2Column);
        cell.setCellStyle(percentageStyle);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Sensitivity LDA:");
        cell = row.createCell(classResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("K"+ldaCorrectRow+"/K"+totalSumRow);
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Sensitivity MS-DIAL");
        cell = row.createCell(adductResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("F"+lbCorrectRow+"/F"+totalSumRow);
        cell = row.createCell(legend2Column);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Sensitivity MS-DIAL:");
        cell = row.createCell(classResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("K"+lbCorrectRow+"/K"+totalSumRow);
        
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Positive Predictive Value LDA:");
        cell = row.createCell(adductResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("F"+ldaCorrectRow+"/(F"+ldaCorrectRow+"+F"+ldaFPRow+")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(percentageStyle);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Positive Predictive Value LDA:");
        cell = row.createCell(classResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("K"+ldaCorrectRow+"/(K"+ldaCorrectRow+"+K"+ldaFPRow+")");

        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Positive Predictive Value MS-DIAL:");
        cell = row.createCell(adductResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("F"+lbCorrectRow+"/(F"+lbCorrectRow+"+F"+lbFPRow+")");
        cell = row.createCell(legend2Column);
        cell.setCellStyle(percentageStyle);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Positive Predictive Value MS-DIAL:");
        cell = row.createCell(classResultColumn);
        cell.setCellStyle(percentageStyle);
        cell.setCellType(Cell.CELL_TYPE_FORMULA);
        cell.setCellFormula("K"+lbCorrectRow+"/(K"+lbCorrectRow+"+K"+lbFPRow+")");

        if (classInfo.getSns()>1){
          rowCount++;
          rowCount++;
          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(leftHeaderStyle);
          cell.setCellValue("Molecular species/structure evaluation");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(leftHeaderStyle);
          cell.setCellValue("Molecular species/structure evaluation - adduct insensitive:");
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          totalSumRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified species:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(H2:H"+finalRowNumberAnalytes+",\"<>*_FP\",H2:H"+finalRowNumberAnalytes+",\"<>\",H2:H"+finalRowNumberAnalytes+",\"<>no structure\",H2:H"+finalRowNumberAnalytes+",\"<>*_noDIAL\")");       
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified species:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(H2:H"+finalRowNumberAnalytes+",\"<>*_FP\",H2:H"+finalRowNumberAnalytes+",\"<>\",H2:H"+finalRowNumberAnalytes+",\"<>no structure\",H2:H"+finalRowNumberAnalytes+",\"<>*_noDIAL\","+ldaMS2Ignore+"2:"+ldaMS2Ignore+finalRowNumberAnalytes+",\"<>TRUE\")");                 
          row =  sheet.createRow(rowCount);
          rowCount++;
          ldaCorrectRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified LDA species:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\">0\")");     
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified LDA species:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\">0\","+ldaMS2Ignore+"2:"+ldaMS2Ignore+finalRowNumberAnalytes+",\"<>TRUE\")");
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          lbCorrectRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified MS-DIAL species:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\">0\")");       
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Total number of identified MS-DIAL species:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\">0\","+dialMS2Ignore+"2:"+dialMS2Ignore+finalRowNumberAnalytes+",\"<>TRUE\")");
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          ldaFPRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False positives LDA:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"<>-1\",K2:K"+finalRowNumberAnalytes+",\"<>0\",K2:K"+finalRowNumberAnalytes+",\"<2\")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False positives LDA:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"<>-1\",K2:K"+finalRowNumberAnalytes+",\"<>0\",K2:K"+finalRowNumberAnalytes+",\"<2\","+ldaMS2Ignore+"2:"+ldaMS2Ignore+finalRowNumberAnalytes+",\"<>TRUE\")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          lbFPRow = rowCount;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False positives MS-DIAL:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"<>-1\",L2:L"+finalRowNumberAnalytes+",\"<>0\",L2:L"+finalRowNumberAnalytes+",\"<2\")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False positives MS-DIAL:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"<>-1\",L2:L"+finalRowNumberAnalytes+",\"<>0\",L2:L"+finalRowNumberAnalytes+",\"<2\","+dialMS2Ignore+"2:"+dialMS2Ignore+finalRowNumberAnalytes+",\"<>TRUE\")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False negatives LDA:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"=-1\")+COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"=-3\")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False negatives LDA:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"=-1\","+ldaMS2Ignore+"2:"+ldaMS2Ignore+finalRowNumberAnalytes+",\"<>TRUE\")+COUNTIFS(K2:K"+finalRowNumberAnalytes+",\"=-3\","+ldaMS2Ignore+"2:"+ldaMS2Ignore+finalRowNumberAnalytes+",\"<>TRUE\")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False negatives MS-DIAL:");
          cell = row.createCell(adductResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"=-1\")+COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"=-3\")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("False negatives MS-DIAL:");
          cell = row.createCell(classResultColumn);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"=-1\","+dialMS2Ignore+"2:"+dialMS2Ignore+finalRowNumberAnalytes+",\"<>TRUE\")+COUNTIFS(L2:L"+finalRowNumberAnalytes+",\"=-3\","+dialMS2Ignore+"2:"+dialMS2Ignore+finalRowNumberAnalytes+",\"<>TRUE\")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Sensitivity LDA:");
          cell = row.createCell(adductResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("F"+ldaCorrectRow+"/F"+totalSumRow);
          cell = row.createCell(legend2Column);
          cell.setCellStyle(percentageStyle);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Sensitivity LDA:");
          cell = row.createCell(classResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("K"+ldaCorrectRow+"/K"+totalSumRow);
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Sensitivity MS-DIAL");
          cell = row.createCell(adductResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("F"+lbCorrectRow+"/F"+totalSumRow);
          cell = row.createCell(legend2Column);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Sensitivity MS-DIAL:");
          cell = row.createCell(classResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("K"+lbCorrectRow+"/K"+totalSumRow);
          
          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Positive Predictive Value LDA:");
          cell = row.createCell(adductResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("F"+ldaCorrectRow+"/(F"+ldaCorrectRow+"+F"+ldaFPRow+")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(percentageStyle);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Positive Predictive Value LDA:");
          cell = row.createCell(classResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("K"+ldaCorrectRow+"/(K"+ldaCorrectRow+"+K"+ldaFPRow+")");

          row =  sheet.createRow(rowCount);
          rowCount++;
          cell = row.createCell(legendColumn);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Positive Predictive Value MS-DIAL:");
          cell = row.createCell(adductResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("F"+lbCorrectRow+"/(F"+lbCorrectRow+"+F"+lbFPRow+")");
          cell = row.createCell(legend2Column);
          cell.setCellStyle(percentageStyle);
          cell.setCellStyle(boldStyle);
          cell.setCellValue("Positive Predictive Value MS-DIAL:");
          cell = row.createCell(classResultColumn);
          cell.setCellStyle(percentageStyle);
          cell.setCellType(Cell.CELL_TYPE_FORMULA);
          cell.setCellFormula("K"+lbCorrectRow+"/(K"+lbCorrectRow+"+K"+lbFPRow+")");
          
        }
        rowCount++;
        rowCount++;
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(leftHeaderStyle);
        cell.setCellValue("Legend:");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(fullyCorrectStyle);
        cell.setCellValue("the name is green if the species is correct");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(notFoundStyle);
        cell.setCellValue("the name is red if it is a false positive");
        rowCount++;
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(fullyCorrectStyle);
        cell.setCellValue("the LDA/MS-DIAL columns are green, if the hit was found");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(ambigStyle);
        cell.setCellValue("the LDA/MS-DIAL column are cyan, if the hit was found additionally at a wrong RT");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(onlyMS1Style);
        cell.setCellValue("the LDA/MS-DIAL column are orange, if the hit was found only at a wrong RT");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(notFoundStyle);
        cell.setCellValue("the LDA/MS-DIAL column are red, if the hit was not reported at all or if the hit is an FP");
        
        rowCount++;
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellStyle(boldStyle);
        cell.setCellValue("Identification Codes:");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("2: correctly identified without any FPs");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("1: correct hit found, but FPs at wrong retention times are reported");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("0: identified only by MS1, or FP not detected");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("-1: false negative hit");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("-2: false positive hit");
        row =  sheet.createRow(rowCount);
        rowCount++;
        cell = row.createCell(legendColumn);
        cell.setCellValue("-3: false negative and false positive hit, e.g., hit at correct RT not found, but something at wrong RT reported");
      }
            
      resultWorkbook.write(out);
      out.close();
    }catch(Exception ex){
      ex.printStackTrace();
    }      
  }

  /**
   *     
   * @param analyte
   * @return int[0] nr of C atoms; int[1] nr of dbs; int[2] nr of OH
   * @throws HydroxylationEncodingException 
   */
  private int[] decodeAnalyte(String analyte) throws HydroxylationEncodingException {
    int[] result = new int[3];
    result[2] = Settings.getLcbHydroxyEncoding().getHydroxyNumber(analyte.substring(0,1));
    analyte = analyte.substring(1);
    result[0] = Integer.parseInt(analyte.substring(0,analyte.indexOf(":")));
    result[1] = Integer.parseInt(analyte.substring(analyte.indexOf(":")+1));
    return result;
  }
  
  
  private LdaLBLASTCompareVO getLdaMSDialSpeciesComparision(String ms1Analyte, LipidClassInfoVO info, LinkedHashMap<String,ReferenceInfoVO> species,
      Vector<LipidParameterSet> ldaAnalytes, Vector<MSDialEntry> ident, Vector<MSDialEntry> dialMS1, LipidomicsAnalyzer analyzer,float mz, String lipidClass,
      boolean useStructure) throws Exception{  
    //split LDA identifications in MS1 identifications and ones where spectra were found
    //this vector contains identifications where no ms2 is present
    Vector<LipidParameterSet> ldaMS1Only = new Vector<LipidParameterSet>();
    //this vector contains all ms2 identifications, including the ones where no structural information can be obtained
    Vector<LipidParameterSet> ldaMs2SpectraIdentified = new Vector<LipidParameterSet>();
    //this hashtable contains identifications where the fatty acids were determined
    Hashtable<String,Vector<LipidomicsMSnSet>> ldaMs2Evidence = new Hashtable<String,Vector<LipidomicsMSnSet>>();
    Hashtable<String,String> ldaOKAccordingToReference = new Hashtable<String,String>();
    Hashtable<String,String> dialOKAccordingToReference = new Hashtable<String,String>();

    boolean isFP = false;
    boolean doNotUseInEvaluation = false;
    for (LipidParameterSet set : ldaAnalytes){
      if (set instanceof LipidomicsMSnSet){
        LipidomicsMSnSet ldaAnalyte = (LipidomicsMSnSet)set;
        ldaMs2SpectraIdentified.add(set);
        boolean ms1Only = ldaAnalyte.getStatus()==2 ? true : false;
        if (!ms1Only) {
          for (Object msnNames : ldaAnalyte.getMSnIdentificationNames()){
            Vector<LipidomicsMSnSet> sameMs2Ident = new Vector<LipidomicsMSnSet>();
            String nameString = null;
            if (msnNames instanceof Vector){
              nameString = ((Vector<String>)msnNames).get(0);
            }else{
              nameString = ((String)msnNames);
            }
            nameString = nameString.replaceAll("/", "_");
            ////if (quant.getAnalyteClass().equalsIgnoreCase("DG")&&set.getModificationName().equalsIgnoreCase("NH4"))nameString+="_-";
            for (String key : ldaMs2Evidence.keySet()){
              if (StaticUtils.isAPermutedVersion(key, nameString,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
                nameString = key;
                sameMs2Ident = ldaMs2Evidence.get(key);
              }
            }
            sameMs2Ident.add(ldaAnalyte);
            ldaMs2Evidence.put(nameString, sameMs2Ident);
          }
        }
      } else ldaMS1Only.add(set);
    }
    
    String correctMS1Name = new String(ms1Analyte);
    int ldaMS1Evidence = 0;
    String ldaMS1Name = "not reported";
    int dialMS1Evidence = 0;
    String dialMS1Name = "not reported";
    String suggestedDialName = "";
    String rtIdeal = "";
    String rtLDA = "";
    String rtDial = "";
    LinkedHashMap<String,String> rtsDial = new LinkedHashMap<String,String>();
    double dialScore = 0d;
    double dialWrongScore = 0d;
    boolean noMS2 = false;
    boolean ms1Only = false;
    Vector<LdaLBLASTCompareVO> ms2Compare = new Vector<LdaLBLASTCompareVO>();
    //now check if LDA or LB identified the analyte at the correct RT
    if (species!=null && species.size()>0){
      doNotUseInEvaluation = true;
      for (ReferenceInfoVO refInfo : species.values()) {
        if (refInfo.useInEvaluation()) {
          doNotUseInEvaluation = false;
          break;
        }
      }
      boolean ms1AnalyteDetected = false;
      for (ReferenceInfoVO ref : species.values()) {
        for (double rt : ref.getCorrectRts()) {
          if (rtIdeal.indexOf(Calculator.FormatNumberToString(rt,2))==-1)
            rtIdeal += Calculator.FormatNumberToString(rt,2)+" ";
        }
      }
      if (rtIdeal.length()>0) rtIdeal = rtIdeal.substring(0,rtIdeal.length()-1);


      // this is the MS1 check of the LDA
      ldaMS1Evidence = -1;
      for (LipidParameterSet set : ldaMs2SpectraIdentified){
        boolean isAtCorrectRt = isHitAtCorrectRt(Double.parseDouble(set.getRt()), info, species);
        if (isAtCorrectRt){
          ms1AnalyteDetected = true;
          if (ldaMS1Evidence==-3 || ldaMS1Evidence==-2) ldaMS1Evidence=1;
          else if (ldaMS1Evidence==-1 || ldaMS1Evidence==0) ldaMS1Evidence=2;
        } else {
          if (ldaMS1Evidence==-2||ldaMS1Evidence==-1) ldaMS1Evidence=-3;
          else if (ldaMS1Evidence==0) ldaMS1Evidence=-2;
          else if (ldaMS1Evidence==2) ldaMS1Evidence=1;          
        }
        rtLDA += set.getRt()+" ";
      }
      if (rtLDA.length()>0) rtLDA = rtLDA.substring(0,rtLDA.length()-1);      
      //this is the ms1 check of MS-DIAL
      dialMS1Evidence = -1;
      String hitRt;
      if (ident!=null){
        for (MSDialEntry detect : ident){
          hitRt = Calculator.FormatNumberToString(detect.getRt(),2d);
          rtsDial.put(hitRt, hitRt);
          boolean isAtCorrectRt = isHitAtCorrectRt(detect.getRt(), info, species);
          if (isAtCorrectRt){
            ms1AnalyteDetected = true;
            if (detect.getScore()>dialScore) {
              dialScore = detect.getScore();
              suggestedDialName = detect.getDialClassName()+" "+detect.getDialMs1Name();
            }
            if (dialMS1Evidence==-3 || dialMS1Evidence==-2) dialMS1Evidence=1;
            else if (dialMS1Evidence==-1 || dialMS1Evidence==0) dialMS1Evidence=2;            
          } else {
            if (dialMS1Evidence==-2||dialMS1Evidence==-1) dialMS1Evidence=-3;
            else if (dialMS1Evidence==0) dialMS1Evidence=-2;
            else if (dialMS1Evidence==2) dialMS1Evidence=1;
            if (dialMS1Evidence==-2 || dialMS1Evidence==-3) {
              if (detect.getScore()>dialWrongScore) {
                dialWrongScore = detect.getScore();
                suggestedDialName = detect.getDialClassName()+" "+detect.getDialMs1Name();
              }              
            }
          }        
        }
      }
      for (String dialRt : rtsDial.keySet()) rtDial += dialRt+" ";
      if (rtDial.length()>0) rtDial = rtDial.substring(0,rtDial.length()-1);
      if ((ldaMS1Evidence==-1||ldaMS1Evidence==-3)  && (dialMS1Evidence==-1||dialMS1Evidence==-3)){
        if (!areThereAnyMSnSpectra(analyzer,mz,info,species)){
          if (ldaMS1Evidence==-3) ldaMS1Evidence=-2;
          if (dialMS1Evidence==-3) dialMS1Evidence=-2;
          if (ldaMS1Evidence==-1) ldaMS1Evidence=0;
          if (dialMS1Evidence==-1) dialMS1Evidence=0;
          noMS2 = true;
        }else{
          ms1AnalyteDetected = true;
        }
      }
      
      //this sets the displayed names for the MS1 identifications for LDA and MS-DIAL
      if (ldaMS1Evidence>0 || ldaMS1Evidence<-1) ldaMS1Name = ms1Analyte;
      else{
        boolean foundMS1Only = false;
        String rt = "";
        for (LipidParameterSet set : ldaMS1Only){
          if (isHitAtCorrectRt(Double.parseDouble(set.getRt()), info, species)){
            rt += set.getRt()+" ";
            foundMS1Only = true;
          }
        }
        if (foundMS1Only){
          ldaMS1Name = ms1Analyte;
          ms1Only = true;
          rt = rt.substring(0,rt.length()-1);
          rtLDA = rt;
          if (dialMS1Evidence<1 && ldaMS1Evidence==-1){
            ms1AnalyteDetected = true;
            ldaMS1Evidence=0;
          }
        }
      }
      if (dialMS1Evidence>0 || dialMS1Evidence<-1) dialMS1Name =  suggestedDialName;//lbMS1Name = ms1Analyte;
      if (!ms1AnalyteDetected && !noMS2){
        correctMS1Name += "_FP";
        isFP = true;
      }
      //here the check for the ms2 evidence starts
      if (info.getSns()>1 && ms1AnalyteDetected && (ldaMS1Evidence>0||dialMS1Evidence>0)){
        for (String key : species.keySet()){
          ReferenceInfoVO ref = species.get(key);
          String keyWOSlash = key.replaceAll("/", "_");
          boolean ms2AnalyteDetected = false;
          String ldaMS2Name = "not reported";
          String rtMS2Ideal = "";
          for (double rt : ref.getCorrectRts()) {
            if (rtMS2Ideal.indexOf(Calculator.FormatNumberToString(rt,2))==-1)
              rtMS2Ideal += Calculator.FormatNumberToString(rt,2)+" ";
          }
          if (rtMS2Ideal.length()>0)
            rtMS2Ideal = rtMS2Ideal.substring(0,rtMS2Ideal.length()-1);
          String rtMS2LDA = "";
          String rtMS2Dial = "";
          // here starts the LDA check of the MS2 evidence
          int ldaMS2Evidence = -1;
          String nameString = null;
          for (String ldaName : ldaMs2Evidence.keySet()){
            Vector<LipidomicsMSnSet> ldaHits = ldaMs2Evidence.get(ldaName);
            if (ldaHits.size()==0) continue;
            String toCompare = new String(ldaName);
            if (lipidClass.equalsIgnoreCase("DG")&&ldaHits.get(0).getModificationName().equalsIgnoreCase("NH4")){
              toCompare += "_-";
            }
                                                                                                                             
            if (!StaticUtils.isAPermutedVersion(keyWOSlash, toCompare,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) continue;
            ldaOKAccordingToReference.put(ldaName, ldaName);
            boolean nameDidNotChange = false;
            for (LipidomicsMSnSet set : ldaHits){              
              // display name for LDA
              for (Object msnNames : set.getMSnIdentificationNames()){
                String oneName = null;
                if (msnNames instanceof Vector){
                  oneName = ((Vector<String>)msnNames).get(0);
                }else{
                  oneName = ((String)msnNames);
                }
                if (!oneName.replaceAll("/","_").equalsIgnoreCase(ldaName))continue;
                if (msnNames instanceof Vector){
                  oneName = "";
                  for (String name: ((Vector<String>)msnNames)){
                    oneName+=name+";";
                  }
                  oneName = oneName.substring(0,oneName.length()-1);
                }else{
                  oneName = ((String)msnNames);
                }
                if (nameString==null){
                  nameDidNotChange = true;
                  nameString = oneName;
                }else if (nameDidNotChange){
                  if (!nameString.equalsIgnoreCase(oneName)){
                    nameDidNotChange = false;
                    nameString = ldaName;
                  }
                }
              }
              
              boolean isAtCorrectRt = isHitAtCorrectRt(Double.parseDouble(set.getRt()), info, ref);
              if (isAtCorrectRt){
                ms2AnalyteDetected = true;
                if (ldaMS2Evidence==-3 || ldaMS2Evidence==-2) ldaMS2Evidence=1;
                else if (ldaMS2Evidence==-1 || ldaMS2Evidence==0) ldaMS2Evidence=2;
              } else {
                if (ldaMS2Evidence==-2||ldaMS2Evidence==-1) ldaMS2Evidence=-3;
                else if (ldaMS2Evidence==0) ldaMS2Evidence=-2;
                else if (ldaMS2Evidence==2) ldaMS2Evidence=1;          
              }
              rtMS2LDA += set.getRt()+" ";
            }
            if (ldaMS2Evidence>0 || ldaMS2Evidence<-1) ldaMS2Name = nameString; 
          }
          if (rtMS2LDA.length()>0) rtMS2LDA.substring(0,rtMS2LDA.length()-1);
          
          // here starts the LB check of the MS2 evidence
          Hashtable<String,String> dialNames = new  Hashtable<String,String>();
          int dialMS2Evidence = -1;
          double dialMs2Score = 0d;
          if (ident!=null && ident.size()>0){
            //for (String lbName : ident.getIdentifications(info.getSns())){
            boolean foundCorrectRt = false;
            for (MSDialEntry entry : ident){
              if (entry.getDialMs2Name()==null || entry.getDialMs2Name().length()==0)
                continue;
              String dialMs2Name = entry.getDialMs2Name();
              if (entry.getLdaMs2Name()!=null && entry.getLdaMs2Name().length()>0)
                dialMs2Name = entry.getLdaMs2Name();
              if (!StaticUtils.isAPermutedVersion(keyWOSlash, dialMs2Name.replaceAll("/", "_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) continue;
              dialNames.put(entry.getDialMs2Name(), entry.getDialMs2Name());
              dialOKAccordingToReference.put(dialMs2Name, dialMs2Name);
              double score = entry.getScore();
              if (!foundCorrectRt && score>dialMs2Score) dialMs2Score=score;
              //for (String rt : ident.getRetentionTimes(lbName, info.getSns())){
              boolean isAtCorrectRt = isHitAtCorrectRt(entry.getRt(), info, ref);
              if (isAtCorrectRt){
                ms2AnalyteDetected = true;
                if (!foundCorrectRt)
                  dialMs2Score=score;
                if (score>dialMs2Score) dialMs2Score=score;                  
                foundCorrectRt=true;
                if (dialMS2Evidence==-3 || dialMS2Evidence==-2) dialMS2Evidence=1;
                else if (dialMS2Evidence==-1 || dialMS2Evidence==0) dialMS2Evidence=2;
              } else {
                if (dialMS2Evidence==-2||dialMS2Evidence==-1) dialMS2Evidence=-3;
                else if (dialMS2Evidence==0) dialMS2Evidence=-2;
                else if (dialMS2Evidence==2) dialMS2Evidence=1;          
              }
              rtMS2Dial += entry.getRt()+" ";
              //}
            }
          }
          String dialName = "";
          if (dialMS2Evidence!=-1 && dialNames.size()>0){
            for (String name : dialNames.keySet()) dialName+=name+" ";
          }
          if (dialName.length()>1)dialName=dialName.substring(0,dialName.length()-1);
          else dialName = "not reported";
          if (rtMS2Dial.length()>1) rtMS2Dial = rtMS2Dial.substring(0,rtMS2Dial.length()-1);
          if (ldaMS2Evidence!=-1 || dialMS2Evidence!=-1){
            if (!useStructure && dialMS2Evidence==-1) {
              dialMS2Evidence = 0;
            }
            ms2Compare.add(new LdaLBLASTCompareVO(key+(!useStructure ? "_noDIAL" : ""),ldaMS2Name,dialName,ldaMS2Evidence,dialMS2Evidence,rtMS2Ideal,rtMS2LDA,rtMS2Dial,false,
                dialMs2Score, mz, false));
          }
        }
      }
    }else{
      //System.out.println("!!!!!!!!!!!! "+ms1Analyte);
      correctMS1Name += "_FP";
      for (LipidParameterSet set : ldaMs2SpectraIdentified){
        ldaMS1Name = ms1Analyte;
        ldaMS1Evidence = LdaLBLASTCompareVO.FALSE_POSITIVE;
      }
      LinkedHashMap<String,String> rts = new LinkedHashMap<String,String>();
      if (ident!=null && ident.size()>0){
        ////dialMS1Name = ms1Analyte;
        dialMS1Evidence = LdaLBLASTCompareVO.FALSE_POSITIVE;
        for (MSDialEntry detect : ident){
          if (detect.getScore()>dialScore) {
            dialScore = detect.getScore();
            dialMS1Name = detect.getDialClassName()+" "+detect.getDialMs1Name();
          }
          String rtString = Calculator.FormatNumberToString(detect.getRt(),2);
          rts.put(rtString,rtString);
        }
      }
      if (rts.size()>0){
        for(String rt : rts.keySet()) rtDial += rt+" ";
        rtDial = rtDial.substring(0,rtDial.length()-1);
      }
      isFP = true;
    }
    //add now the ms2 FPs of the LDA including the ones that MS-DIAL shows too
    if (info.getSns()>1 && (ldaMS1Evidence<-1 || ldaMS1Evidence>0)){
      for (String ldaKey : ldaMs2Evidence.keySet()){
        if (ldaOKAccordingToReference.containsKey(ldaKey)) continue;
        Vector<LipidomicsMSnSet> ldaHits = ldaMs2Evidence.get(ldaKey);
        String nameString = null;
        if (ldaHits.size()==0) continue;
        boolean nameDidNotChange = false;
        String rtMS2LDA = "";
        for (LipidomicsMSnSet set : ldaHits){
          for (Object msnNames : set.getMSnIdentificationNames()){
            String oneName = null;
            if (msnNames instanceof Vector){
              oneName = ((Vector<String>)msnNames).get(0);
            }else{
              oneName = ((String)msnNames);
            }
            if (!oneName.replaceAll("/","_").equalsIgnoreCase(ldaKey))continue;
            if (msnNames instanceof Vector){
              oneName = "";
              for (String name: ((Vector<String>)msnNames)){
                oneName+=name+";";
              }
              oneName = oneName.substring(0,oneName.length()-1);
            }else{
              oneName = ((String)msnNames);
            }
            if (nameString==null){
              nameDidNotChange = true;
              nameString = oneName;
            }else if (nameDidNotChange){
              if (!nameString.equalsIgnoreCase(oneName)){
                nameDidNotChange = false;
                nameString = ldaKey;
              }
            }
            rtMS2LDA += set.getRt()+" ";
          }
        }
        int dialMS2Evidence = 0;
        double dialMs2Prob = 0d;
        String rtMS2Dial = "";
        Hashtable<String,String> dialNames = new  Hashtable<String,String>();
        //check if MS-DIAL finds the FP too
        if (ident!=null && (dialMS1Evidence<-1  || dialMS1Evidence>0)){
          //for (String lbNm : ident.getIdentifications(info.getSns())){
          for (MSDialEntry entry : ident){
            if (entry.getDialMs2Name()==null || entry.getDialMs2Name().length()==0)
              continue;
            String dialMs2Name = entry.getDialMs2Name();
            if (entry.getLdaMs2Name()!=null && entry.getLdaMs2Name().length()>0)
              dialMs2Name = entry.getLdaMs2Name();
            if (!StaticUtils.isAPermutedVersion(ldaKey, dialMs2Name.replaceAll("/", "_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)) continue;
            dialMS2Evidence = LdaLBLASTCompareVO.FALSE_POSITIVE;
            dialNames.put(entry.getDialMs2Name(), entry.getDialMs2Name());
            dialOKAccordingToReference.put(dialMs2Name, dialMs2Name);
            double score = entry.getScore();
            if (score>dialMs2Prob) dialMs2Prob=score;
            //for (String rt : ident.getRetentionTimes(lbNm, info.getSns())){
            String rtString = Calculator.FormatNumberToString(entry.getRt(),2);
            rtMS2Dial += rtString+" ";              
            //}
          }         
        }
        String dialName = "";
        if (dialMS2Evidence!=-1 && dialNames.size()>0){
          for (String name : dialNames.keySet()) dialName+=name+" ";
        }
        if (dialName.length()>1)dialName=dialName.substring(0,dialName.length()-1);
        if (rtMS2Dial.length()>1) rtMS2Dial = rtMS2Dial.substring(0,rtMS2Dial.length()-1);
        ms2Compare.add(new LdaLBLASTCompareVO(nameString+"_FP",nameString,dialName,-2,dialMS2Evidence,"",rtMS2LDA,rtMS2Dial,true,
            dialMs2Prob, mz, false));
      }
    }
    //add now the ms2 FPs that were found by MS-DIAL only
    if (info.getSns()>1 && ident!=null && (dialMS1Evidence<-1  || dialMS1Evidence>0)){
      int dialMS2Evidence = 0;
      Hashtable<String,String> dialNames = new  Hashtable<String,String>();
      Hashtable<String, Hashtable<String,String>> dialNameVariations = new Hashtable<String, Hashtable<String,String>>();
      Hashtable<String,String> rtsMS2Dial = new Hashtable<String,String>();
      Hashtable<String,Double> scores = new Hashtable<String,Double>();
      //for (String lbNm : ident.getIdentifications(info.getSns())){
      for (MSDialEntry entry : ident){
        if (entry.getDialMs2Name()==null || entry.getDialMs2Name().length()==0)
          continue;
        String dialMs2Name = entry.getDialMs2Name();
        if (entry.getLdaMs2Name()!=null && entry.getLdaMs2Name().length()>0)
          dialMs2Name = entry.getLdaMs2Name();
        if (dialOKAccordingToReference.containsKey(dialMs2Name)) continue;
        dialMS2Evidence = LdaLBLASTCompareVO.FALSE_POSITIVE;
        Hashtable<String,String> nameVariations = new Hashtable<String,String>();
        String rts = "";
        double score = 0d;
        String key = dialMs2Name;
        for (String other : dialNames.keySet()){
          if (StaticUtils.isAPermutedVersion(other.replaceAll("/", "_"),dialMs2Name.replaceAll("/","_"),LipidomicsConstants.CHAIN_SEPARATOR_NO_POS)){
            nameVariations = dialNameVariations.get(other);
            
            rts = rtsMS2Dial.get(other);
            key = other;
            break;
          }
        }
        nameVariations.put(entry.getDialMs2Name(), entry.getDialMs2Name());
        double scr = entry.getScore();
        if (scr>score) score = scr;
        //for (String rt : ident.getRetentionTimes(lbNm, info.getSns())){
        String rtString = Calculator.FormatNumberToString(entry.getRt(),2);
        rts += rtString+" ";              
        //}
        dialNames.put(key, dialMs2Name);
        dialNameVariations.put(key, nameVariations);
        rtsMS2Dial.put(key, rts);
        scores.put(key, score);
      }
      if (dialMS2Evidence==LdaLBLASTCompareVO.FALSE_POSITIVE && dialNames.size()>0){
        for (String key : dialNames.keySet()){
          String dialName = "";
          Hashtable<String,String> nameVariations = dialNameVariations.get(key);
          String rtMS2Dial = rtsMS2Dial.get(key);
          double dialMs2Score = scores.get(key);
          if (dialMS2Evidence!=-1 && dialNames.size()>0){
            for (String name : nameVariations.keySet()) dialName+=name+" ";
          }
          if (dialName.length()>1)dialName=dialName.substring(0,dialName.length()-1);
          if (rtMS2Dial.length()>0) rtMS2Dial = rtMS2Dial.substring(0,rtMS2Dial.length()-1);
          ms2Compare.add(new LdaLBLASTCompareVO(dialName+"_FP","not reported",dialName,0,dialMS2Evidence,"","",rtMS2Dial,true,
              dialMs2Score, mz, false));
        }
      }
    }
    LdaLBLASTCompareVO compare = null;
    boolean createCompareVO = true;
    if (ldaMS1Evidence==0 && dialMS1Evidence==0 && !ms1Only) createCompareVO = false;
    if (ldaMS1Evidence<1 && dialMS1Evidence<1 && noMS2) {
      correctMS1Name += "_noMS2";
      if (dialMS1Evidence==0 && dialMS1!=null && dialMS1.size()>0) {
        String hitRt;
        for (MSDialEntry detect : dialMS1) {
          hitRt = Calculator.FormatNumberToString(detect.getRt(),2d);
          rtsDial.put(hitRt, hitRt);
          boolean isAtCorrectRt = isHitAtCorrectRt(detect.getRt(), info, species);
          if (isAtCorrectRt){
            dialMS1Name = detect.getDialClassName()+" "+detect.getDialMs1Name();
            break;
          }        
        }
      }
    }
    if (createCompareVO){
      if (doNotUseInEvaluation)
        dialMS1Evidence = 0;
      compare = new LdaLBLASTCompareVO(correctMS1Name,ldaMS1Name,dialMS1Name,ldaMS1Evidence,dialMS1Evidence,rtIdeal,rtLDA,rtDial,isFP,
        dialScore, mz, ms1Only);
      if (ms2Compare.size()>0) compare.setMs2Evidence(ms2Compare);
        
    }
    return compare;
  }
}
