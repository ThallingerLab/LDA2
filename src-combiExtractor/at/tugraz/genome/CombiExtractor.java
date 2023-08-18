package at.tugraz.genome;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.Hashtable;
import java.util.Vector;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.util.CellRangeAddress;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.Calculator;

public class CombiExtractor
{

  //specification of column headers shown
  private final static String COLUMN_FILE = "file";
  //removed, since I should simply take here the lipid class
  //private final static String COLUMN_SPECIES = "TG";
  private final static String COLUMN_RETENTION_TIME = "RT";
  private final static String COLUMN_FORMULA = "formula";
  private final static String COLUMN_MASS = "m/z";
  private final static String COLUMN_MOLSPECIES = "Chains Composition";
  private final static String COLUMN_LDA_RESULT = "LDA [%]";

  
  public static void main(String[] args)
  {
    new CombiExtractor();  
  }
  
  public CombiExtractor() {
    convertLDAResultCombisToOverviewExel();
  }
  
  private void convertLDAResultCombisToOverviewExel() {
    
    //specification of column sequence
    int fileColumn = 0;
    int speciesColumn = 1;
    int rtColumn = 2;
    int formulaColumn = 3;
    int mzColumn = 4;
    int molSpeciesColumn = 5;
    int firstChainColumn = 6;
    
    File dir = new File(System.getProperty("user.dir"));
    String resultFilename = "CombinationOverview.xlsx";
    System.out.println(dir);
    String outFilename = dir+File.separator+resultFilename;
    BufferedOutputStream out = null;
    Workbook resultWorkbook = null;
    try {
      out = new BufferedOutputStream(new FileOutputStream(outFilename));
      resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      
      String longestLipidClass = "";
      String longestChainFragment = "";
      int longestFile = 0;
      int longestSpecies = 0;
      int longestRt = 0;
      int longestFormula = 0;
      String longestAdduct = "";
      int longestMass = 0;
      int longestMolSpecies = 0;
      int longestChain = 0;
      int longestFragmentValue = 0;
      int mostLeftColumn = 0;
      
      File[] files = dir.listFiles();
      for (File file : files){
        if (!file.getName().endsWith(".xlsx") || file.getName().equalsIgnoreCase(resultFilename)) continue;
        QuantificationResult returnParam = LDAResultReader.readResultFile(file.getName(), new Hashtable<String,Boolean>());
        if (returnParam==null){
          System.out.println("It was not possible to translate the file: "+file.getName());
          continue;
        }
        String sheetName = file.getName();
        if (sheetName.length()>30)
          sheetName = sheetName.substring(0, 30);
        Sheet sheet = resultWorkbook.createSheet(sheetName);        
        
        //TODO: header and explanatory row shall be created for each lipid class
        int rowCount = 0;
        
        for (String lipidClass : returnParam.getIdentifications().keySet()){
          Vector<LipidParameterSet> params = returnParam.getIdentifications().get(lipidClass);
          @SuppressWarnings("rawtypes")
          Vector chainAndFrag = extractChainAndFragmentInfo(params);
          @SuppressWarnings("unchecked")
          Vector<String> adducts = (Vector<String>)chainAndFrag.get(0);
          @SuppressWarnings("unchecked")
          Hashtable<String,Integer> nrOfChains = (Hashtable<String,Integer>)chainAndFrag.get(1);
          @SuppressWarnings("unchecked")
          Hashtable<String,Hashtable<String,Double>> strongestAdduct = (Hashtable<String,Hashtable<String,Double>>)chainAndFrag.get(2);
          for (String adduct : adducts) {
            int numberOfChains = nrOfChains.get(adduct);
            if (numberOfChains<2)
              continue;
            String strongestChainFragment = strongestAdduct.get(adduct).keySet().iterator().next();
            if (strongestChainFragment.length()>longestChainFragment.length())
              longestChainFragment = strongestChainFragment;
          
            //explanatory row
            Row row = sheet.createRow(rowCount);
            createCell(row, headerStyle, firstChainColumn, "m/z");
            sheet.addMergedRegion(new CellRangeAddress(rowCount, rowCount, firstChainColumn,firstChainColumn+numberOfChains-1));
            rowCount++;
          
            //header row
            row = sheet.createRow(rowCount);
            createCell(row, headerStyle, fileColumn, COLUMN_FILE);
            createCell(row, headerStyle, speciesColumn, lipidClass);
            if (lipidClass.length()>longestLipidClass.length())
              longestLipidClass = lipidClass;
            createCell(row, headerStyle, rtColumn, COLUMN_RETENTION_TIME);
            createCell(row, headerStyle, formulaColumn, COLUMN_FORMULA);
            String massAndAdduct = COLUMN_MASS+"_"+adduct;
            createCell(row, headerStyle, mzColumn, massAndAdduct);
            if (massAndAdduct.length()>longestAdduct.length())
              longestAdduct = massAndAdduct;
            createCell(row, headerStyle, molSpeciesColumn, COLUMN_MOLSPECIES);
            for (int i=firstChainColumn; i!=(firstChainColumn+numberOfChains); i++) {
              createCell(row, headerStyle, i, strongestChainFragment);
            }
            createCell(row, headerStyle, firstChainColumn+numberOfChains, COLUMN_LDA_RESULT);
            if (firstChainColumn+numberOfChains>mostLeftColumn)
              mostLeftColumn = firstChainColumn+numberOfChains;
            rowCount++;

            for (LipidParameterSet set:params) {
              if (!set.getModificationName().equalsIgnoreCase(adduct) || !(set instanceof LipidomicsMSnSet))
                continue;
              LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
              if (msn.getStatus()<LipidomicsMSnSet.FRAGMENTS_DETECTED)
                continue;
              //now everything is fine -> print out the information
              row = sheet.createRow(rowCount);
              createCell(row, null, fileColumn, file.getName());
              if (file.getName().length()>longestFile)
                longestFile = file.getName().length();
              createCell(row, null, speciesColumn, msn.getNameStringWithoutRt());
              if (msn.getNameStringWithoutRt().length()>longestSpecies)
                longestSpecies = msn.getNameStringWithoutRt().length();
              createCell(row, null, rtColumn, msn.getRt());
              if (msn.getRt().length()>longestRt)
                longestRt = msn.getRt().length();
              Hashtable<String,Integer> categorized = StaticUtils.categorizeFormula(msn.getAnalyteFormula());
              String formula = StaticUtils.getFormulaInHillNotation(categorized, false);
              createCell(row, null, formulaColumn, formula);
              if (formula.length()>longestFormula)
                longestFormula = formula.length();
              String mzString = Calculator.FormatNumberToString((double)msn.Mz[0],3);
              createCell(row, null, mzColumn, mzString);
              if (mzString.length()>longestMass)
                longestMass = mzString.length();
              Vector<String> sortedCombis = new Vector<String>();
              for (String combi : msn.getChainCombinationRelativeAreas().keySet()) {
                double relArea = msn.getChainCombinationRelativeAreas().get(combi);
                boolean wasHigher=false;
                for (int i=0; i!=sortedCombis.size(); i++) {
                  if (relArea>msn.getChainCombinationRelativeAreas().get(sortedCombis.get(i))) {
                    wasHigher=true;
                    sortedCombis.add(i, combi);
                    break;
                  }
                }
                if (!wasHigher)
                  sortedCombis.add(combi);
              }
              boolean isFirst = true;
              for (String combi : sortedCombis) {
                if (!isFirst)
                  row = sheet.createRow(rowCount);
                else
                  isFirst = false;
                String sortedCombiName = "";
                Vector<FattyAcidVO> sorted = new Vector<FattyAcidVO>();
                try {
                	Vector<FattyAcidVO> decoded = StaticUtils.decodeLipidNamesFromChainCombi(combi);
                	sorted = StaticUtils.sortChainVOs(decoded);
                	sortedCombiName = StaticUtils.getHumanReadableCombiName(sorted, returnParam.getFaHydroxyEncoding(), returnParam.getLcbHydroxyEncoding());
                } catch (LipidCombinameEncodingException ex) {
                	ex.printStackTrace();
                }
                createCell(row, null, molSpeciesColumn, sortedCombiName);
                if (sortedCombiName.length()>longestMolSpecies)
                  longestMolSpecies = sortedCombiName.length();
                int count = 0;
                for (FattyAcidVO fa : sorted) {
                	count++;
                	String encoded = StaticUtils.encodeLipidNameForCreatingCombis(fa, true);
                	Hashtable<String,CgProbe> fragments = msn.getChainFragments().get(encoded);
                  if (!fragments.containsKey(strongestChainFragment))
                    continue;
                  String mz = Calculator.FormatNumberToString((double)fragments.get(strongestChainFragment).Mz,3);
                  createCell(row, null, molSpeciesColumn+count, mz);
                  if (mz.length()>longestFragmentValue)
                    longestFragmentValue = mz.length();
                }
                String percentValue = String.valueOf(Calculator.roundDBL(msn.getChainCombinationRelativeAreas().get(combi)*100d,2,BigDecimal.ROUND_HALF_UP));
                createCell(row, null, molSpeciesColumn+count+1, percentValue);
                if (percentValue.length()>longestFragmentValue)
                  longestFragmentValue = percentValue.length();               
                rowCount++;
              }
            }                   
          }
        }
        
      }
      for (int i=0; i!=resultWorkbook.getNumberOfSheets(); i++) {
        Sheet sheet = resultWorkbook.getSheetAt(i);
        setColumnWidth(sheet, fileColumn, COLUMN_FILE, longestFile);
        setColumnWidth(sheet, speciesColumn, longestLipidClass, longestSpecies);
        setColumnWidth(sheet, rtColumn, COLUMN_RETENTION_TIME, longestRt);
        setColumnWidth(sheet, formulaColumn, COLUMN_FORMULA, longestFormula);
        setColumnWidth(sheet, mzColumn, longestAdduct, longestMass);
        setColumnWidth(sheet, molSpeciesColumn, COLUMN_MOLSPECIES, longestMolSpecies);

        for (int j=firstChainColumn; j!=(mostLeftColumn+1); j++) {
          setColumnWidth(sheet, j, longestChainFragment, longestFragmentValue);
        }

      }
      
      resultWorkbook.write(out);
      resultWorkbook.close();
      out.close();
    } catch (IOException | ExcelInputFileException | ChemicalFormulaException ex) {
      ex.printStackTrace();
    } finally {
      if (resultWorkbook!=null) {
        try {resultWorkbook.close();}catch (IOException e) {e.printStackTrace();}
      }
      if (out!=null) {
        try {out.close();}catch (IOException e) { e.printStackTrace();}
      }
    }
  }
 
  
  
  protected static CellStyle getHeaderStyle(Workbook wb){
    CellStyle arial12style = wb.createCellStyle();
    Font arial12font = wb.createFont();
    arial12font.setBoldweight(Font.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(CellStyle.ALIGN_CENTER);
    return arial12style;
  }

  @SuppressWarnings({ "rawtypes", "unchecked" })
  private Vector extractChainAndFragmentInfo(Vector<LipidParameterSet> params) {
    Vector result = new Vector();
    Hashtable<String,Integer> nrOfChains = new Hashtable<String,Integer>();
    Hashtable<String,Hashtable<String,Double>> strongestAdduct = new Hashtable<String,Hashtable<String,Double>>();
    //int nrOfChains = 0;
    int nr;
    Hashtable<String,Double> strongestFragmentInfo;
    for (LipidParameterSet set : params) {
      if (!(set instanceof LipidomicsMSnSet))
        continue;
      LipidomicsMSnSet msn = (LipidomicsMSnSet)set;
      if (msn.getStatus()<LipidomicsMSnSet.FRAGMENTS_DETECTED)
        continue;
      String adduct = msn.getModificationName();
      String strongestFragment;
      double highestArea;
      if (nrOfChains.containsKey(adduct)) {
        nr = nrOfChains.get(adduct);
        strongestFragment = strongestAdduct.get(adduct).keySet().iterator().next();
        highestArea = strongestAdduct.get(adduct).get(strongestFragment);
      }else{
        nr = 0;
        strongestFragment = null;
        highestArea = 0;
      }
      try {
      	for (Object obj : msn.getMSnIdentificationNames()){
          if (obj instanceof Vector){
            if (getNumberOfChainsFromName((String)((Vector)obj).get(0))>nr)
              nr = getNumberOfChainsFromName((String)((Vector)obj).get(0));
          }else{
            if (getNumberOfChainsFromName((String)obj)>nr)          
              nr = getNumberOfChainsFromName((String)obj);
          }
        }
      } catch (LipidCombinameEncodingException ex) {
      	ex.printStackTrace();
      }
      for (String chain : msn.getChainFragments().keySet()) {
        for (String fragment : msn.getChainFragments().get(chain).keySet()) {
          if (msn.getChainFragments().get(chain).get(fragment).Area>highestArea) {
            strongestFragment = fragment;
            highestArea = msn.getChainFragments().get(chain).get(fragment).Area;
          }
        }
      }
      nrOfChains.put(adduct, nr);
      strongestFragmentInfo = new Hashtable<String,Double>();
      strongestFragmentInfo.put(strongestFragment, highestArea);
      strongestAdduct.put(adduct,  strongestFragmentInfo);
    }
    Vector<String> adductsSorted = new Vector<String>();
    boolean added = false;
    String otherAdduct;
    for (String adduct : strongestAdduct.keySet()) {
      for (int i=0; i!=adductsSorted.size(); i++) {
        otherAdduct = adductsSorted.get(i);
        if (strongestAdduct.get(adduct).values().iterator().next()>strongestAdduct.get(otherAdduct).values().iterator().next()) {
          adductsSorted.add(i,adduct);
          added = true;
          break;
        }
      }
      if (!added)
        adductsSorted.add(adduct);
    }
    
    result.add(adductsSorted);
    result.add(nrOfChains);
    result.add(strongestAdduct);    
    return result;
  }
  
  private int getNumberOfChainsFromName(String name) {
    String[] splitted = name.replaceAll("/", "_").split("_");
    int nrOfChains = splitted.length;
    for (String chain : splitted) {
      if (chain.equalsIgnoreCase("-"))
        nrOfChains--;
    }
    return nrOfChains;
  }
  
  private Cell createCell(Row row, CellStyle style, int pos, String value){
    Cell cell = row.createCell(pos);
    cell.setCellValue(value);
    if (cell!=null) cell.setCellStyle(style);
    return cell;
  }

  private void setColumnWidth(Sheet sheet, int column, String headerValue, int longestValue){
    int columnWidth = (int)((headerValue.length()*256)*ExcelUtils.BOLD_MULT);
    if ((longestValue+1)*256>columnWidth) columnWidth =  (longestValue+1)*256;
    sheet.setColumnWidth(column,columnWidth); 
  }

}
