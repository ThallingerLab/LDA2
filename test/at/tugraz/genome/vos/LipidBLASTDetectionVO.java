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

package at.tugraz.genome.vos;

import java.util.Hashtable;
import java.util.StringTokenizer;

import at.tugraz.genome.exception.LipidBLASTException;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidBLASTDetectionVO
{
  private String lipidClass_;
  private String adduct_;
  private String ms1Name_;
  private String ms2Name_;
  private String sn1_;
  private String sn2_;
  private String sn3_;
  private String retentionTime_;
  private int scanNum_;
  private String fileName_;
  private String precursorMass_;
  private int rank_;
  private String library_;
  private long libraryId_;
  private String mass_;
  private String deltaMz_;
  private String libMz_;
  private int score_;
  private int dotProduct_;
  private double probability_;
  private int revDot_;
  
  private static Hashtable<String,String> classLookup_ = new Hashtable<String,String>();
  private static Hashtable<String,String> adductLookup_ = new Hashtable<String,String>();
  
  public LipidBLASTDetectionVO(String ms1Name, String adduct,
      String ms2Name, String retentionTime, int scanNum,
      String fileName, String precursorMass, int rank, String library,
      long libraryId, String mass, String deltaMz, String libMz, int score,
      int dotProduct, double probability, int revDot) throws LipidBLASTException
  {
    initLookups();
    String[] classAndName = getClassAndName(ms1Name);
    this.lipidClass_ = classAndName[0];
    this.adduct_ = adduct;
    if (adductLookup_.containsKey(adduct_)) this.adduct_ = adductLookup_.get(adduct_);
    // error of LipidBlast - it should be [M+H]+ but is [M]+
    else if (lipidClass_.equalsIgnoreCase("SM") && adduct_.equalsIgnoreCase("[M]+")) this.adduct_ = "H";
    this.ms1Name_ = classAndName[1];
    if (lipidClass_.equalsIgnoreCase("Cer")){
      this.ms1Name_=ms2Name.substring(ms2Name.indexOf("/")+1,ms2Name.lastIndexOf(")"));
      if (this.ms1Name_.indexOf("(")!=-1) this.ms1Name_=this.ms1Name_.substring(0,this.ms1Name_.indexOf("("));
    } else if (lipidClass_.equalsIgnoreCase("SM")){
      int cAtoms = Integer.parseInt(ms1Name_.substring(0,ms1Name_.indexOf(":")))-18;
      int dbs = Integer.parseInt(ms1Name_.substring(ms1Name_.indexOf(":")+1))-1;
      if (dbs<0) dbs=100;
      this.ms1Name_ = String.valueOf(cAtoms)+":"+dbs;
    }
    this.ms2Name_ = ms2Name;
    if (lipidClass_.equalsIgnoreCase("Cer")){
      this.ms2Name_ = "Cer ("+ms1Name_+"/0:0)";
    } else if (lipidClass_.equalsIgnoreCase("SM")){
      this.ms2Name_ = "SM ("+ms1Name_+")";
    }
    this.assignMSnPositions();
    this.retentionTime_ = retentionTime;
    this.scanNum_ = scanNum;
    this.fileName_ = fileName;
    this.precursorMass_ = precursorMass;
    this.rank_ = rank;
    this.library_ = library;
    this.libraryId_ = libraryId;
    this.mass_ = mass;
    this.deltaMz_ =  deltaMz;
    this.libMz_ = libMz;
    this.score_ = score;
    this.dotProduct_ = dotProduct;
    this.probability_ = probability;
    this.revDot_ = revDot;
  }
  
  private void initLookups(){
    classLookup_.put("plasmenyl-PC", "P-PC");
    classLookup_.put("plasmenyl-PE", "P-PE");
    classLookup_.put("lysoPC", "LPC");
    classLookup_.put("lysoPE", "LPE");
    adductLookup_.put("[M+H]+", "H");
    adductLookup_.put("[M]+", "H");
    adductLookup_.put("[M+Na]+", "Na");
    adductLookup_.put("[M+NH4]+", "NH4");
    adductLookup_.put("[M-H]-","-H");
    adductLookup_.put("[M+HCOO]-","HCOO");
    adductLookup_.put("[M-CH3]-","-CH3");
  }
  
  private String[] getClassAndName(String ms1Name) throws LipidBLASTException{
    String[] classAndName = new String[2];
    StringTokenizer tokenizer = new StringTokenizer(ms1Name," ");
    if (tokenizer.countTokens()==2){
      classAndName[0] = tokenizer.nextToken();
      classAndName[1] = tokenizer.nextToken();
    } else if (tokenizer.countTokens()==1 && ms1Name.startsWith("N-") && ms1Name.endsWith("-sphing-4-enine")){
      classAndName[0] = "Cer";
    } else throw new LipidBLASTException("A lipid class name must have only one empty space in the name; the name \""+ms1Name+"\" has not!");
    if (classLookup_.containsKey(classAndName[0])) classAndName[0] = classLookup_.get(classAndName[0]);
    return classAndName;
  }
  
  private void assignMSnPositions(){
    String fas = this.ms2Name_.substring(this.ms2Name_.indexOf("(")+1,this.ms2Name_.lastIndexOf(")"));
    StringTokenizer tokenizer = new  StringTokenizer(fas,"/");
    int tokenNumber = 0;
    sn1_ = "-";
    sn2_ = "-";
    sn3_ = "-";
    while (tokenizer.hasMoreTokens()){
      tokenNumber++;
      String fa = tokenizer.nextToken();
      if (fa.indexOf("(")!=-1) fa = fa.substring(0,fa.indexOf("(")).trim();
      if (fa.equalsIgnoreCase("0:0")) fa = "-";
      if (tokenNumber==1) sn1_ = fa;
      if (tokenNumber==2) sn2_ = fa;
      if (tokenNumber==3) sn3_ = fa;
    }
  }

  public String getLipidClass()
  {
    return lipidClass_;
  }

  public String getAdduct()
  {
    return adduct_;
  }

  public String getMs1Name()
  {
    return ms1Name_;
  }

  public String getMs2Name()
  {
    return ms2Name_;
  }

  public String getSn1()
  {
    return sn1_;
  }

  public String getSn2()
  {
    return sn2_;
  }

  public String getSn3()
  {
    return sn3_;
  }

  public String getRetentionTime()
  {
    return retentionTime_;
  }

  public int getScanNum()
  {
    return scanNum_;
  }

  public String getFileName()
  {
    return fileName_;
  }

  public String getPrecursorMass()
  {
    return precursorMass_;
  }

  public int getRank()
  {
    return rank_;
  }

  public String getLibrary()
  {
    return library_;
  }

  public long getLibraryId()
  {
    return libraryId_;
  }

  public String getMass()
  {
    return mass_;
  }

  public String getDeltaMz()
  {
    return deltaMz_;
  }

  public String getLibMz()
  {
    return libMz_;
  }

  public int getScore()
  {
    return score_;
  }

  public int getDotProduct()
  {
    return dotProduct_;
  }

  public double getProbability()
  {
    return probability_;
  }

  public int getRevDot()
  {
    return revDot_;
  }
  
  public String toString(){
    String result = "";
    result += "Class: "+lipidClass_+"; ";
    result += "Adduct: "+adduct_+"; ";
    result += "MS1: "+ms1Name_+"; ";
    result += "MS2: "+sn1_+"/"+sn2_+"/"+sn3_+"; ";
    result += "RT: "+retentionTime_+"; ";
    result += "Rank: "+rank_+"; ";
    result += "Prob: "+probability_+"; ";
    result += "Scan: "+scanNum_+"; ";
    result += "Prec: "+precursorMass_+"; ";
    result += "Mass: "+mass_+"; ";
    result += "LibMz: "+libMz_+"; ";
    result += "Delta: "+deltaMz_+"; ";
    result += "Lib: "+library_+"; ";
    result += "Lib-ID: "+libraryId_+"; ";
    result += "Score: "+score_+"; ";
    result += "Dot: "+dotProduct_+"; ";
    result += "Rev-Dot: "+revDot_+"; ";
    result += "File: "+fileName_+"; ";
    return result;
  }
  
  
}
