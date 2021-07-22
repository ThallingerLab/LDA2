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

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.LDAResultReader;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.FloatStringVO;
import at.tugraz.genome.voutils.GeneralComparator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LDAToFASummaryConverter extends ConverterBase
{

  protected final static int COLUMN_FILENAME = 0;
  
  protected final static String TEXT_FILENAME = "File name";
  
  private String dirName_;
  
  public LDAToFASummaryConverter(String dirName){
    dirName_ = dirName;
  }
  
  public void convert() throws ExcelInputFileException, FileNotFoundException, LipidCombinameEncodingException{
    File dir = new File(dirName_);
    File[] files = dir.listFiles();
    
    //reading the results and adding them to a hash map
    //and creating a list of available lipid classes
    LinkedHashMap<String,QuantificationResult> rawResults = new LinkedHashMap<String,QuantificationResult>();
    LinkedHashSet<String> lipidClasses = new LinkedHashSet<String>();
    for (File file : files){
      String fileName = file.getName();
      if (!fileName.endsWith(".xlsx")) continue;
      
      QuantificationResult returnParam = LDAResultReader.readResultFile(file.getAbsolutePath(), new Hashtable<String,Boolean>());
      if (returnParam==null){
        System.out.println("It was not possible to translate the file: "+fileName);
        return;
      }
      Hashtable<String,Vector<LipidParameterSet>> results = returnParam.getIdentifications();
      if (results == null || results.keySet().size()==0)
        continue;
      boolean containsData = false;
      for (Vector<LipidParameterSet> idents : results.values()){
        if (idents.size()>0){
          containsData = true;
          break;
        }
      }
      if (!containsData){
        System.out.println("The file does not contain any data: "+fileName);
        continue;
      }
      rawResults.put(fileName, returnParam);
      for (String lipidClass : results.keySet()){
        if (!lipidClasses.contains(lipidClass))
          lipidClasses.add(lipidClass);
      }
    }
    
    writeResults(lipidClasses,rawResults);
  }
  
  private void writeResults(LinkedHashSet<String> lipidClasses, LinkedHashMap<String,QuantificationResult> rawResults)
      throws FileNotFoundException, LipidCombinameEncodingException{
    String outFilename = dirName_+File.separator+"Summary.xlsx";
    BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outFilename));
    Vector<HydroxyEncoding> ohEncodings = unifyHydroxyEncoding(rawResults);
    try{
      Workbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      for (String lClass : lipidClasses){
        Sheet sheet = resultWorkbook.createSheet(lClass);
        LongestCount lCount = new LongestCount();
        Row row = sheet.createRow(lCount.rowCount_);
        lCount.rowCount_++;
        
        this.createCell(row, headerStyle, COLUMN_FILENAME, TEXT_FILENAME);
        this.createCell(row, headerStyle, COLUMN_LIPID_SPECIES+1, TEXT_LIPID_SPECIES);
        this.createCell(row, headerStyle, COLUMN_ADDUCT+1, TEXT_ADDUCT);
        this.createCell(row, headerStyle, COLUMN_RT+1, TEXT_RT);
        this.createCell(row, headerStyle, COLUMN_LIPID_SPECIES_INTENSITY+1, TEXT_LIPID_SPECIES_INTENSITY);
        this.createCell(row, headerStyle, COLUMN_CHAIN+1, TEXT_CHAIN);
        this.createCell(row, headerStyle, COLUMN_CHAIN_PERCENT+1, TEXT_CHAIN_PERCENT);
        this.createCell(row, headerStyle, COLUMN_CHAIN_INTENSITY+1, TEXT_CHAIN_INTENSITY);
        this.createCell(row, headerStyle, COLUMN_MOLECULAR_SPECIES+1, TEXT_MOLECULAR_SPECIES);

        Vector<Vector<String>> bothInfo = getSortedSpeciesListAndAdducts(lClass, rawResults);
        Vector<String> sortedSpecies = bothInfo.get(0);
        Vector<String> mods = bothInfo.get(1);
        String species;
        String faReadable;
        Vector<FattyAcidVO> chains;
        for (String speciesEncoded : sortedSpecies){
          chains = StaticUtils.decodeLipidNamesFromChainCombi(speciesEncoded);
          species = StaticUtils.getHumanReadableChainName(chains.get(0), ohEncodings.get(0), ohEncodings.get(1), StaticUtils.areThereOhInCombi(chains));
          //first add species without fatty acids;
          writeExcelRowsOfOneSpecies(lCount,sheet,species,null,rawResults.keySet(),mods,getAreasWithoutChainInfo(lClass,species,rawResults));
          //second add species with fatty acids;
          LinkedHashMap<String,Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>>> sortedFAs = getSortedFattyAcids(lClass,species,speciesEncoded,rawResults);
          for (String fa: sortedFAs.keySet()) {
            chains = StaticUtils.decodeLipidNamesFromChainCombi(fa);
            faReadable = StaticUtils.getHumanReadableChainName(chains.get(0), ohEncodings.get(0), ohEncodings.get(1), StaticUtils.areThereOhInCombi(chains));
            writeExcelRowsOfOneSpecies(lCount,sheet,species,faReadable,rawResults.keySet(),mods,sortedFAs.get(fa));
          }
        }
        setColumnWidth(sheet, COLUMN_FILENAME, TEXT_FILENAME, lCount.longestFilename_);
        setColumnWidth(sheet, COLUMN_LIPID_SPECIES+1, TEXT_LIPID_SPECIES, lCount.longestSpecies_);
        setColumnWidth(sheet, COLUMN_ADDUCT+1, TEXT_ADDUCT, lCount.longestAdduct_);
        setColumnWidth(sheet, COLUMN_RT+1, TEXT_RT, lCount.longestRt_);
        setColumnWidth(sheet, COLUMN_LIPID_SPECIES_INTENSITY+1, TEXT_LIPID_SPECIES_INTENSITY, lCount.longestSpeciesIntensity_);
        setColumnWidth(sheet, COLUMN_CHAIN+1, TEXT_CHAIN, lCount.longestChain_);
        setColumnWidth(sheet, COLUMN_CHAIN_PERCENT+1, TEXT_CHAIN_PERCENT, lCount.longestChainPercent_);
        setColumnWidth(sheet, COLUMN_CHAIN_INTENSITY+1, TEXT_CHAIN_INTENSITY, lCount.longestChainIntensity_);
        setColumnWidth(sheet, COLUMN_MOLECULAR_SPECIES+1, TEXT_MOLECULAR_SPECIES, lCount.longestMolecularSpecies_);
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
  
  @SuppressWarnings("unchecked")
  private Vector<Vector<String>> getSortedSpeciesListAndAdducts(String lClass, LinkedHashMap<String,QuantificationResult> rawResults) throws LipidCombinameEncodingException{
    Vector<Vector<String>> bothInfo = new Vector<Vector<String>>();
    Set<String> species = new HashSet<String>();
    Hashtable<String,FloatStringVO> modAreas = new Hashtable<String,FloatStringVO>();
    for (QuantificationResult res : rawResults.values()){
      if (!res.getIdentifications().containsKey(lClass))
        continue;
      for (LipidParameterSet set : res.getIdentifications().get(lClass)){
        species.add(StaticUtils.decodeHumanReadableChain(set.getNameStringWithoutRt(), res.getFaHydroxyEncoding(), res.getLcbHydroxyEncoding(), res.getConstants().isAlexTargetlist()).getChainId());
        if (!modAreas.containsKey(set.getModificationName()))
          modAreas.put(set.getModificationName(), new FloatStringVO(set.getModificationName(),0f));
        modAreas.get(set.getModificationName()).addValue(set.getArea());
      }
    }
    bothInfo.add(sortFAsInAscendingOrder(species));
    List<FloatStringVO> mods = new ArrayList<FloatStringVO>(modAreas.values());
    Collections.sort(mods,new GeneralComparator("at.tugraz.genome.lda.vos.FloatStringVO", "getValue", "java.lang.Float"));
    Vector<String> modsSorted = new Vector<String>();
    for (int i=(mods.size()-1);i!=-1;i--){
      modsSorted.add(mods.get(i).getKey());
    }
    bothInfo.add(modsSorted);
    return bothInfo;
  }
  
  
  private LinkedHashMap<String,Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>>> getSortedFattyAcids(String lClass, String speciesName, String speciesEncoded, LinkedHashMap<String,QuantificationResult> rawResults) throws LipidCombinameEncodingException{
    //first key: chain sorted; second key: modification; third key: experiment, fourth key: retention time; value: details
    Hashtable<String,Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>>> details = new Hashtable<String,Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>>>();
    for (String exp : rawResults.keySet()){
      QuantificationResult res  = rawResults.get(exp);
      if (!res.getIdentifications().containsKey(lClass))
        continue;
      for (LipidParameterSet set : res.getIdentifications().get(lClass)){
        if (!set.getNameStringWithoutRt().equalsIgnoreCase(speciesName))
          continue;
        if (!StaticUtils.isThereChainInformationAvailable(set.getNameStringWithoutRt(), set))
          continue;
        float ms1Area = getMS1Area(set);
        String ms1AreaString = String.valueOf(ms1Area);
        LipidomicsMSnSet msn = (LipidomicsMSnSet) set;
        LinkedHashMap<String,String[]> faDetails = extractFaIntensityDetails(msn, ms1AreaString, ms1Area);
        for (String fa : faDetails.keySet()){
          if (!details.containsKey(fa))
            details.put(fa, new Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>>());
          if (!details.get(fa).containsKey(msn.getModificationName()))
            details.get(fa).put(msn.getModificationName(), new Hashtable<String,LinkedHashMap<String,String[]>>());
          if (!details.get(fa).get(msn.getModificationName()).containsKey(exp))
            details.get(fa).get(msn.getModificationName()).put(exp, new LinkedHashMap<String,String[]>());
          details.get(fa).get(msn.getModificationName()).get(exp).put(msn.getRt(), faDetails.get(fa));          
        }
      }
    }
    LinkedHashMap<String,Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>>> detailsSorted = new LinkedHashMap<String,Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>>>();
    for (String fa : sortFAsInAscendingOrder(details.keySet()))
      detailsSorted.put(fa, details.get(fa));
    return detailsSorted;
  }
  
  private Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>> getAreasWithoutChainInfo(String lClass, String speciesName, LinkedHashMap<String,QuantificationResult> rawResults) throws LipidCombinameEncodingException{
    //first key: modification; second key: experiment, third key: retention time; value: details
    Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>> areasWithoutChainInfo = new Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>>();
    for (String exp : rawResults.keySet()){
      QuantificationResult res  = rawResults.get(exp);
      if (!res.getIdentifications().containsKey(lClass))
        continue;
      for (LipidParameterSet set : res.getIdentifications().get(lClass)){
        if (!set.getNameStringWithoutRt().equalsIgnoreCase(speciesName))
          continue;
        if (StaticUtils.isThereChainInformationAvailable(set.getNameStringWithoutRt(), set))
          continue;
        float ms1Area = getMS1Area(set);
        String ms1AreaString = String.valueOf(ms1Area);
        String[] details = new String[1];
        details[0] = ms1AreaString;
        if (!areasWithoutChainInfo.containsKey(set.getModificationName()))
          areasWithoutChainInfo.put(set.getModificationName(), new Hashtable<String,LinkedHashMap<String,String[]>>());
        if (!areasWithoutChainInfo.get(set.getModificationName()).containsKey(exp))
          areasWithoutChainInfo.get(set.getModificationName()).put(exp, new LinkedHashMap<String,String[]>());
        areasWithoutChainInfo.get(set.getModificationName()).get(exp).put(set.getRt(), details);
      }
    }
    return areasWithoutChainInfo;
  }
  
  private void writeExcelRowsOfOneSpecies(LongestCount lCount, Sheet sheet, String species, String chain, Set<String> expSequence, Vector<String> mods,
      Hashtable<String,Hashtable<String,LinkedHashMap<String,String[]>>> details){
    Row row;
    for (String mod :  mods){
      if (!details.containsKey(mod))
        continue;
      for (String exp : expSequence){
        if (!details.get(mod).containsKey(exp))
          continue;
        LinkedHashMap<String,String[]> detailsOfExp =  details.get(mod).get(exp);
        for (String rt : detailsOfExp.keySet()){
          String[] detailsOfOneHit = detailsOfExp.get(rt);
          row = sheet.createRow(lCount.rowCount_);
          lCount.rowCount_++;
          createCell(row, null, COLUMN_FILENAME, exp);
          if (exp.length()>lCount.longestFilename_) lCount.longestFilename_ = exp.length();
          createCell(row, null, COLUMN_LIPID_SPECIES+1, species);
          if (species.length()>lCount.longestSpecies_) lCount.longestSpecies_ = species.length();
          createCell(row, null, COLUMN_ADDUCT+1, mod);
          if (mod.length()>lCount.longestAdduct_) lCount.longestAdduct_ = mod.length();
          createNumericCell(row, null, COLUMN_RT+1, rt);
          if (rt.length()>lCount.longestRt_) lCount.longestRt_ = rt.length();
          createNumericCell(row, null, COLUMN_LIPID_SPECIES_INTENSITY+1, detailsOfOneHit[0]);
          if (detailsOfOneHit[0].length()>lCount.longestSpeciesIntensity_) lCount.longestSpeciesIntensity_ = detailsOfOneHit[0].length();
          if (chain!=null){
            createCell(row, null, COLUMN_CHAIN+1, chain);
            if (chain.length()>lCount.longestChain_) lCount.longestChain_ = chain.length();
            createNumericCell(row, null, COLUMN_CHAIN_PERCENT+1, detailsOfOneHit[1]);
            if (detailsOfOneHit[1].length()>lCount.longestChainPercent_) lCount.longestChainPercent_ = detailsOfOneHit[1].length();
            createNumericCell(row, null, COLUMN_CHAIN_INTENSITY+1, detailsOfOneHit[2]);
            if (detailsOfOneHit[2].length()>lCount.longestChainIntensity_) lCount.longestChainIntensity_ = detailsOfOneHit[2].length();
            createCell(row, null, COLUMN_MOLECULAR_SPECIES+1, detailsOfOneHit[3]);
            if (detailsOfOneHit[3].length()>lCount.longestMolecularSpecies_) lCount.longestMolecularSpecies_ = detailsOfOneHit[3].length();

          }
        }
      }
    }
  }
  
  private class LongestCount {
    
    public int rowCount_ = 0;
    public int longestFilename_ = 0;
    public int longestSpecies_ = 0;
    public int longestAdduct_ = 0;
    public int longestRt_ = 0;
    public int longestSpeciesIntensity_ = 0;
    public int longestChain_ = 0;
    public int longestChainPercent_ = 0;
    public int longestChainIntensity_ = 0;
    public int longestMolecularSpecies_ = 0;
    
    public LongestCount(){
      
    }
  }
  
  private Vector<HydroxyEncoding> unifyHydroxyEncoding(LinkedHashMap<String,QuantificationResult> rawResults){
    Vector<HydroxyEncoding> encodings = new Vector<HydroxyEncoding>();
    for (QuantificationResult res : rawResults.values()) {
      if (encodings.size()==0) {
        encodings.add(res.getFaHydroxyEncoding());
        encodings.add(res.getLcbHydroxyEncoding());
      } else {
        encodings.get(0).mergeHydroxyEncodings(res.getFaHydroxyEncoding());
        encodings.get(1).mergeHydroxyEncodings(res.getLcbHydroxyEncoding());
      }
    }
    return encodings;
  }
/*  
  results[0] = ms1AreaString;
  results[1] = percent;
  results[2] = faAreaString;
  results[3] = combiName;*/

}
