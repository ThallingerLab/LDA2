/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2018 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

import java.util.Collection;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Vector;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.utils.ExcelUtils;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */

public abstract class ConverterBase
{
  
  protected final static int COLUMN_LIPID_SPECIES = 0;
  protected final static int COLUMN_ADDUCT = 1;
  protected final static int COLUMN_RT = 2;
  protected final static int COLUMN_LIPID_SPECIES_INTENSITY = 3;
  protected final static int COLUMN_CHAIN = 4;
  protected final static int COLUMN_CHAIN_PERCENT = 5;
  protected final static int COLUMN_CHAIN_INTENSITY = 6;
  protected final static int COLUMN_MOLECULAR_SPECIES = 7;
  
  protected final static String TEXT_LIPID_SPECIES = "Lipid species";
  protected final static String TEXT_ADDUCT = "Adduct";
  protected final static String TEXT_RT = "RT";
  protected final static String TEXT_LIPID_SPECIES_INTENSITY = "Intensity";
  protected final static String TEXT_CHAIN = "Chain";
  protected final static String TEXT_CHAIN_PERCENT = "Percent";
  protected final static String TEXT_CHAIN_INTENSITY = "Intensity";
  protected final static String TEXT_MOLECULAR_SPECIES = "Molecular species";
  
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
  
  protected void setColumnWidth(Sheet sheet, int column, String headerValue, int longestValue){
    int columnWidth = (int)((headerValue.length()*256)*ExcelUtils.BOLD_MULT);
    if ((longestValue+1)*256>columnWidth) columnWidth =  (longestValue+1)*256;
    sheet.setColumnWidth(column,columnWidth); 
  }
  
  protected Cell createNumericCell(Row row, CellStyle style, int pos, String value){
    Cell cell = createNumericCell(row, style, pos, Double.valueOf(value));
    return cell;
  }
  
  protected Cell createNumericCell(Row row, CellStyle style, int pos, double value){
    Cell cell = createCell(row, style, pos, null);
    cell.setCellType(Cell.CELL_TYPE_NUMERIC);
    cell.setCellValue(value);
    return cell;
  }
  
  protected Cell createCell(Row row, CellStyle style, int pos, String value){
    Cell cell = row.createCell(pos);
    cell.setCellValue(value);
    if (cell!=null) cell.setCellStyle(style);
    return cell;
  }

  protected Vector<String> sortFAsInAscendingOrder(Collection<String> unsorted) throws LipidCombinameEncodingException{
    Vector<FattyAcidVO> fas = new Vector<FattyAcidVO>();
    for (String name : unsorted)
      fas.add(StaticUtils.decodeLipidNamesFromChainCombi(name).get(0));
    fas = StaticUtils.sortChainVOs(fas);
    Vector<String> sorted = new Vector<String>();
    for (FattyAcidVO fa : fas) 
      sorted.add(fa.getChainId());
    return sorted;
  }
  
  protected boolean isFaPresent(String fa, String combiName) throws LipidCombinameEncodingException{
   Vector<String> fas = StaticUtils.splitChainCombiToEncodedStrings(combiName,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    for (String oneFa : fas){
      if (oneFa.equalsIgnoreCase(fa))
        return true;
    }
    return false;
  }

  @SuppressWarnings("unchecked")
  private String[] extractFaIntensityDetails(String fa, Vector<String> sortedFAs, LipidomicsMSnSet msn,
      String ms1AreaString, float ms1Area, float faIntensity, float totalIntensity) throws LipidCombinameEncodingException{
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
        if (isFaPresent(fa, (String)combi)) toAdd = getFASortedCombiName((String)combi);
      }else if (combi instanceof Vector){
        for (String onePoss : (Vector<String>)combi){
          if (!isFaPresent(fa, (String)onePoss)) continue;
          if (toAdd.length()>0) toAdd+="|";
          toAdd += getFASortedCombiName(onePoss);
        }
      }
      if (toAdd.length()>0){
        if (combiNames.length()!=0) combiNames.append(";");
        combiNames.append(toAdd);
      }
    }
    String combiName = combiNames.toString();

    String[] results = new String[4];
    results[0] = ms1AreaString;
    results[1] = percent;
    results[2] = faAreaString;
    results[3] = combiName;
    return results;
  }
  
  private String getFASortedCombiName(String name) throws LipidCombinameEncodingException{
    if (name.indexOf("/")==-1 && name.indexOf(LipidomicsConstants.CHAIN_COMBI_SEPARATOR)!=-1)
      return StaticUtils.sortFASequenceUnassigned(name,LipidomicsConstants.CHAIN_SEPARATOR_NO_POS);
    return name;
  }
  
  protected LinkedHashMap<String,String[]> extractFaIntensityDetails(LipidomicsMSnSet msn, String ms1AreaString, float ms1Area) throws LipidCombinameEncodingException{
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
    LinkedHashMap<String,String[]> details = new LinkedHashMap<String,String[]>();
    Vector<String> sortedFAs = sortFAsInAscendingOrder(faIntensities.keySet());
    for (String fa : sortedFAs){
      String[] areaResults = extractFaIntensityDetails(fa, sortedFAs, msn,
          ms1AreaString, ms1Area, faIntensities.get(fa), totalIntensity);
      details.put(fa,areaResults);
    }            
    return details;  
  }
  
  protected float getMS1Area(LipidParameterSet set){
    return set.getArea(0);
  }
  
}
