package at.tugraz.genome;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.LDAResultReader;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.Calculator;

public class LDAToFaConverter
{
  
  private String fileName_;
  
  public LDAToFaConverter(String fileName){
    fileName_ = fileName;
  }
  
  private final static int COLUMN_LIPID_SPECIES = 0;
  private final static int COLUMN_ADDUCT = 1;
  private final static int COLUMN_RT = 2;
  private final static int COLUMN_LIPID_SPECIES_INTENSITY = 3;
  private final static int COLUMN_CHAIN = 4;
  private final static int COLUMN_CHAIN_PERCENT = 5;
  private final static int COLUMN_CHAIN_INTENSITY = 6;
  private final static int COLUMN_MOLECULAR_SPECIES = 7;
  
  private final static String TEXT_LIPID_SPECIES = "Lipid species";
  private final static String TEXT_ADDUCT = "Adduct";
  private final static String TEXT_RT = "RT";
  private final static String TEXT_LIPID_SPECIES_INTENSITY = "Intensity";
  private final static String TEXT_CHAIN = "Chain";
  private final static String TEXT_CHAIN_PERCENT = "Percent";
  private final static String TEXT_CHAIN_INTENSITY = "Intensity";
  private final static String TEXT_MOLECULAR_SPECIES = "Molecular species";
  
  @SuppressWarnings("unchecked")
  public void convert() throws ExcelInputFileException, FileNotFoundException{
    QuantificationResult returnParam = LDAResultReader.readResultFile(fileName_, new Hashtable<String,Boolean>());
    if (returnParam==null){
      System.out.println("It was not possible to translate the file: "+fileName_);
      return;
    }
    Hashtable<String,Vector<LipidParameterSet>> results = returnParam.getIdentifications();
    if (results == null || results.keySet().size()==0){
      System.out.println("The file does not contain any data: "+fileName_);
      return;      
    }
    boolean allZero = true;
    for (Vector<LipidParameterSet> sets : results.values()){
      if (sets.size()!=0){
        allZero = false;
        break;
      }
    }
    if (allZero){
      System.out.println("The file does not contain any data: "+fileName_);
      return;      
    }
    String outFilename = fileName_.substring(0,fileName_.lastIndexOf("."))+"_FA"+fileName_.substring(fileName_.lastIndexOf("."));
    BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outFilename));
    try{
      Workbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      for (String className : results.keySet()){
        Sheet sheet = resultWorkbook.createSheet(className);
        int rowCount = 0;
        Row row = sheet.createRow(rowCount);
        rowCount++;

        int longestSpecies = 0;
        int longestAdduct = 0;
        int longestRt = 0;
        int longestSpeciesIntensity = 0;
        int longestChain = 0;
        int longestChainPercent = 0;
        int longestChainIntensity = 0;
        int longestMolecularSpecies = 0;
        
        this.createCell(row, headerStyle, COLUMN_LIPID_SPECIES, TEXT_LIPID_SPECIES);
        this.createCell(row, headerStyle, COLUMN_ADDUCT, TEXT_ADDUCT);
        this.createCell(row, headerStyle, COLUMN_RT, TEXT_RT);
        this.createCell(row, headerStyle, COLUMN_LIPID_SPECIES_INTENSITY, TEXT_LIPID_SPECIES_INTENSITY);
        this.createCell(row, headerStyle, COLUMN_CHAIN, TEXT_CHAIN);
        this.createCell(row, headerStyle, COLUMN_CHAIN_PERCENT, TEXT_CHAIN_PERCENT);
        this.createCell(row, headerStyle, COLUMN_CHAIN_INTENSITY, TEXT_CHAIN_INTENSITY);
        this.createCell(row, headerStyle, COLUMN_MOLECULAR_SPECIES, TEXT_MOLECULAR_SPECIES);

        
        for (LipidParameterSet set : results.get(className)){
          String displayName = set.getNameStringWithoutRt();
          String rt = set.getRt();
          float ms1Area = 0f;
          // intensity of MS1
          for (CgProbe probe : set.getIsotopicProbes().get(0))
            ms1Area += probe.Area;
          String ms1AreaString = String.valueOf(ms1Area);
          String mod = String.valueOf(set.getModificationName());
          if (set instanceof LipidomicsMSnSet && ((LipidomicsMSnSet)set).getStatus()>=LipidomicsMSnSet.FRAGMENTS_DETECTED){
            LipidomicsMSnSet msn = (LipidomicsMSnSet) set;
            Hashtable<String,Hashtable<String,CgProbe>> chainFrags = msn.getChainFragments();
            float totalIntensity = 0f;
            Hashtable<String,Float> faIntensities = new Hashtable<String,Float>();
            for (String fa : chainFrags.keySet()){
              Hashtable<String,CgProbe> fragments = chainFrags.get(fa);
              float faArea = 0f;
              for (CgProbe fragment : fragments.values()){
                faArea += fragment.Area;
              }
              totalIntensity += faArea;
              faIntensities.put(fa, faArea);
            }            
            Vector<String> sortedFAs = sortFAsInAscendingOrder(faIntensities.keySet());
            for (String fa : sortedFAs){
              float faIntensity = faIntensities.get(fa);
              String percent = "100.00";
              String faAreaString = ms1AreaString;
              if (sortedFAs.size()!=1){
                float relativeValue = faIntensity/totalIntensity;
                percent = String.valueOf(Calculator.roundFloat(100f*relativeValue, 2));
                faAreaString = String.valueOf(ms1Area*relativeValue);
              }
              StringBuilder combiNames = new StringBuilder();
              for (Object combi : msn.getMSnIdentificationNames()){
                String toAdd = "";
                if (combi instanceof String){
                  if (isFaPresent(fa, (String)combi)) toAdd = (String)combi;
                }else if (combi instanceof Vector){
                  for (String onePoss : (Vector<String>)combi){
                    if (!isFaPresent(fa, (String)onePoss)) continue;
                    if (toAdd.length()>0) toAdd+="|";
                    toAdd += onePoss;
                  }
                }
                if (toAdd.length()>0){
                  if (combiNames.length()!=0) combiNames.append(";");
                  combiNames.append(toAdd);
                }
              }
              String combiName = combiNames.toString();
              
              row = sheet.createRow(rowCount);
              rowCount++;
              createCell(row, null, COLUMN_LIPID_SPECIES, displayName);
              if (displayName.length()>longestSpecies) longestSpecies = displayName.length();
              createCell(row, null, COLUMN_ADDUCT, mod);
              if (mod.length()>longestAdduct) longestAdduct = mod.length();
              createNumericCell(row, null, COLUMN_RT, rt);
              if (rt.length()>longestRt) longestRt = rt.length();
              createNumericCell(row, null, COLUMN_LIPID_SPECIES_INTENSITY, ms1AreaString);
              if (ms1AreaString.length()>longestSpeciesIntensity) longestSpeciesIntensity = ms1AreaString.length();
              createCell(row, null, COLUMN_CHAIN, fa);
              if (fa.length()>longestChain) longestChain = fa.length();
              createNumericCell(row, null, COLUMN_CHAIN_PERCENT, percent);
              if (percent.length()>longestChainPercent) longestChainPercent = percent.length();
              createNumericCell(row, null, COLUMN_CHAIN_INTENSITY, faAreaString);
              if (faAreaString.length()>longestChainIntensity) longestChainIntensity = faAreaString.length();
              createCell(row, null, COLUMN_MOLECULAR_SPECIES, combiName);
              if (combiName.length()>longestMolecularSpecies) longestMolecularSpecies = combiName.length();
            }            
          } else {
            row = sheet.createRow(rowCount);
            rowCount++;
            createCell(row, null, COLUMN_LIPID_SPECIES, displayName);
            if (displayName.length()>longestSpecies) longestSpecies = displayName.length();
            createCell(row, null, COLUMN_ADDUCT, mod);
            if (mod.length()>longestAdduct) longestAdduct = mod.length();
            createNumericCell(row, null, COLUMN_RT, rt);
            if (rt.length()>longestRt) longestRt = rt.length();
            createNumericCell(row, null, COLUMN_LIPID_SPECIES_INTENSITY, ms1AreaString);
            if (ms1AreaString.length()>longestSpeciesIntensity) longestSpeciesIntensity = ms1AreaString.length();

          }
        }
        
        setColumnWidth(sheet, COLUMN_LIPID_SPECIES, TEXT_LIPID_SPECIES, longestSpecies);
        setColumnWidth(sheet, COLUMN_ADDUCT, TEXT_ADDUCT, longestAdduct);
        setColumnWidth(sheet, COLUMN_RT, TEXT_RT, longestRt);
        setColumnWidth(sheet, COLUMN_LIPID_SPECIES_INTENSITY, TEXT_LIPID_SPECIES_INTENSITY, longestSpeciesIntensity);
        setColumnWidth(sheet, COLUMN_CHAIN, TEXT_CHAIN, longestChain);
        setColumnWidth(sheet, COLUMN_CHAIN_PERCENT, TEXT_CHAIN_PERCENT, longestChainPercent);
        setColumnWidth(sheet, COLUMN_CHAIN_INTENSITY, TEXT_CHAIN_INTENSITY, longestChainIntensity);
        setColumnWidth(sheet, COLUMN_MOLECULAR_SPECIES, TEXT_MOLECULAR_SPECIES, longestMolecularSpecies);
      }
      resultWorkbook.write(out);
      resultWorkbook.close();
      out.close();
    } catch (IOException iox) {
      try {
        out.close();
      } catch (IOException e) {
        e.printStackTrace();
      }      
    } finally{
      try {
        out.close();
      }
      catch (IOException e) {
        e.printStackTrace();
      }
    }
  }
  
  private static CellStyle getHeaderStyle(Workbook wb){
    CellStyle arial12style = wb.createCellStyle();
    Font arial12font = wb.createFont();
    arial12font.setBoldweight(Font.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(CellStyle.ALIGN_CENTER);
    return arial12style;
  }
  
  private void setColumnWidth(Sheet sheet, int column, String headerValue, int longestValue){
    int columnWidth = (int)((headerValue.length()*256)*ExcelUtils.BOLD_MULT);
    if ((longestValue+1)*256>columnWidth) columnWidth =  (longestValue+1)*256;
    sheet.setColumnWidth(column,columnWidth); 
  }
  
  private Cell createNumericCell(Row row, CellStyle style, int pos, String value){
    Cell cell = createNumericCell(row, style, pos, Double.valueOf(value));
    return cell;
  }
  
  private Cell createNumericCell(Row row, CellStyle style, int pos, double value){
    Cell cell = createCell(row, style, pos, null);
    cell.setCellType(Cell.CELL_TYPE_NUMERIC);
    cell.setCellValue(value);
    return cell;
  }
  
  private Cell createCell(Row row, CellStyle style, int pos, String value){
    Cell cell = row.createCell(pos);
    cell.setCellValue(value);
    if (cell!=null) cell.setCellStyle(style);
    return cell;
  }

  private Vector<String> sortFAsInAscendingOrder(Collection<String> unsorted){
    Vector<String> sorted = new Vector<String>();
    List<Integer> carbonNumbers = new ArrayList<Integer>();
    Hashtable<Integer,Hashtable<Integer,List<String>>> carbonHash = new Hashtable<Integer,Hashtable<Integer,List<String>>>();
    for (String fa : unsorted){
      int prefix = 0;
      char[] chars = fa.toCharArray();
      while (!Character.isDigit(chars[prefix]))
        prefix++;
      int carbons = Integer.parseInt(fa.substring(prefix,fa.indexOf(":")));
      int doubleBonds = Integer.parseInt(fa.substring(fa.indexOf(":")+1));
      Hashtable<Integer,List<String>> sameCarbons = new Hashtable<Integer,List<String>>();
      if (carbonHash.containsKey(carbons)){
        sameCarbons = carbonHash.get(carbons);
      } else {
        carbonNumbers.add(carbons);        
      }
      List<String> sameCandDb = new ArrayList<String>();
      if (sameCarbons.containsKey(doubleBonds)) sameCandDb = sameCarbons.get(doubleBonds);
      sameCandDb.add(fa);
      sameCarbons.put(doubleBonds, sameCandDb);
      carbonHash.put(carbons, sameCarbons);
    }
    Collections.sort(carbonNumbers);
    for (Integer carbon : carbonNumbers){
      Hashtable<Integer,List<String>> sameCarbons = carbonHash.get(carbon);
      List<Integer> doubleBonds = new ArrayList<Integer>(sameCarbons.keySet());
      Collections.sort(doubleBonds);
      for (Integer db : doubleBonds){
        List<String> sameCandDb = sameCarbons.get(db);
        Collections.sort(sameCandDb);
        for (String fa : sameCandDb){
          sorted.add(fa);
        }
      }
    }
    return sorted;
  }
  
  private boolean isFaPresent(String fa, String combiName){
    String[] fas = LipidomicsMSnSet.getFAsFromCombiName(combiName);
    for (String oneFa : fas){
      if (oneFa.equalsIgnoreCase(fa))
        return true;
    }
    return false;
  }
}
