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
package at.tugraz.genome.vos;

import java.util.Hashtable;

import at.tugraz.genome.exception.MSDialException;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MSDialEntry
{
  
  private static Hashtable<String,String> adductLookup_ = new Hashtable<String,String>(){{
    put("[M-H]-","-H");
    put("[M+HCOO]-","HCOO");
    put("[M+H]+","H");
    put("[M+Na]+","Na");
    put("[M+H-H2O]+","-OH");
  }};

  private String id_;
  private String name_;
  private String dialClassName_;
  private String dialMs1Name_;
  private String dialMs2Name_;
  private String ldaClassName_;
  private String ldaMs1Name_;
  private String ldaMs2Name_;
  private int scans_;
  private float rtStart_;
  private float rt_;
  private float rtStop_;
  private double mz_;
  private float area_;
  private String adduct_;
  private String isotope_;
  private float score_;
  private float signalNoise_;
  private String msmsSpectrum_;
  
  
  public MSDialEntry(String id, String name, int scans, float rtStart, float rt, float rtStop, double mz,
      float area, String adduct, String isotope, float score, float signalNoise, String msmsSpectrum) throws MSDialException
  {
    super();
    this.id_ = id;
    this.name_ = name;
    this.categorizeName(mz,rt);
    this.scans_ = scans;
    this.rtStart_ = rtStart;
    this.rt_ = rt;
    this.rtStop_ = rtStop;
    this.mz_ = mz;
    this.area_ = area;
    if (adductLookup_.containsKey(adduct))
      this.adduct_ = adductLookup_.get(adduct);
    else
      throw new MSDialException("The adduct \""+adduct+"\" cannot be found in the lookup table at MSDialEntry.java - please add it!");
    this.isotope_ = isotope;
    this.score_ = score;
    this.signalNoise_ = signalNoise;
    this.msmsSpectrum_ = msmsSpectrum;
  }
  
  
  private void categorizeName(double mz, float rt) throws MSDialException{
    dialClassName_ = null;
    dialMs2Name_ = null;
    ldaClassName_ = null;
    ldaMs1Name_ = null;
    ldaMs2Name_ = null;
    dialMs1Name_ = name_;
    if (dialMs1Name_.startsWith("w/o MS2:"))
      dialMs1Name_ = dialMs1Name_.substring("w/o MS2:".length());
//    System.out.println("name: "+dialMs1Name_);
    if (dialMs1Name_.indexOf(";")!=-1) {
      dialMs2Name_ = dialMs1Name_.substring(dialMs1Name_.indexOf(";")+1);
      dialMs1Name_ = dialMs1Name_.substring(0,dialMs1Name_.indexOf(";"));
    }
    if (dialMs1Name_.indexOf(" ")==-1)
      throw new MSDialException("There is no class in the identification: "+name_);    
    dialClassName_ = dialMs1Name_.substring(0,dialMs1Name_.indexOf(" "));
    dialMs1Name_ = dialMs1Name_.substring(dialClassName_.length()).trim();
    if (dialMs2Name_!=null && dialMs2Name_.length()>0)
      dialMs2Name_ = dialMs2Name_.substring(dialClassName_.length()+1).trim();
    
    ldaClassName_ = dialClassName_;
    ldaMs1Name_ = dialMs1Name_;
    ldaMs2Name_ = dialMs2Name_;
    boolean faHydroxylated = false;
    if (dialClassName_.lastIndexOf("-")!=-1) {
      ldaClassName_ = dialClassName_.substring(0,dialClassName_.lastIndexOf("-"));
      String hydroxClass = dialClassName_.substring(dialClassName_.lastIndexOf("-")+1);
      if (hydroxClass.startsWith("A")||hydroxClass.startsWith("B")||hydroxClass.startsWith("EO")) {
        faHydroxylated = true;
        hydroxClass = hydroxClass.substring(1);
        if (hydroxClass.startsWith("EO"))
          hydroxClass = hydroxClass.substring(1);
      }else if (hydroxClass.startsWith("N")){
        hydroxClass = hydroxClass.substring(1);
      }else if (hydroxClass.startsWith("H")){
        hydroxClass = hydroxClass.substring(1);
        System.out.println("I do not know what the hydroxyClass \""+hydroxClass+"\" might be: "+name_+"   "+mz+";"+rt);        
      }else
        throw new MSDialException("The hydroxylation encoding cannot be decoded: "+name_);
      int nrOfHydroxies = 0;
      if (hydroxClass.equalsIgnoreCase("DS")||hydroxClass.equalsIgnoreCase("S")) {
        nrOfHydroxies = 2;
      }else if (hydroxClass.equalsIgnoreCase("P")) {
        nrOfHydroxies = 3;
      }
      if (faHydroxylated)
        nrOfHydroxies++;
      String hEncoding = "";
      if (nrOfHydroxies==1)
        hEncoding = "m";
      else if (nrOfHydroxies==2)
        hEncoding = "d";
      else if (nrOfHydroxies==3)
        hEncoding = "t";
      else if (nrOfHydroxies==4)
        hEncoding = "q";
      else
        throw new MSDialException("This amount of hydroxylation encodings is not possible in MS-DIAL: "+nrOfHydroxies+"  "+name_);
      ldaMs1Name_ = hEncoding+dialMs1Name_.substring(1);
      if (ldaMs1Name_.endsWith("+O"))
        ldaMs1Name_ = ldaMs1Name_.substring(0,ldaMs1Name_.length()-2);
      else if (ldaMs1Name_.endsWith("+1O"))
        ldaMs1Name_ = ldaMs1Name_.substring(0,ldaMs1Name_.length()-3);
      
    } 
    if (dialMs2Name_!=null){
      ldaMs2Name_ = dialMs2Name_.substring(0,dialMs2Name_.indexOf("/")+1);
      if (faHydroxylated)
        ldaMs2Name_ += "h";
      else
        ldaMs2Name_ += "n";
      ldaMs2Name_ += dialMs2Name_.substring(dialMs2Name_.indexOf("/")+1);
      if (ldaMs2Name_.endsWith("+O"))
        ldaMs2Name_ = ldaMs2Name_.substring(0,ldaMs2Name_.length()-2);
    }
    if (ldaClassName_.equalsIgnoreCase("CerP")) {
      ldaClassName_ = "Cer1P";
      ldaMs1Name_ = "d"+ldaMs1Name_;
    } else if (dialClassName_.equalsIgnoreCase("Sphinganine")||dialClassName_.equalsIgnoreCase("Sphingosine")) {
      ldaClassName_ = "SphBase";
      ldaMs1Name_ = "d"+ldaMs1Name_;
    } else if (dialClassName_.equalsIgnoreCase("Phytosphingosine")) {
      ldaClassName_ = "SphBase";
      ldaMs1Name_ = "t"+ldaMs1Name_;
    }

//    System.out.println(dialClassName_+" ! "+dialMs1Name_+" ! "+dialMs2Name_);
//    System.out.println(ldaClassName_+" ! "+ldaMs1Name_+" ! "+ldaMs2Name_);
    
  }


  public String getId()
  {
    return id_;
  }


  public String getName()
  {
    return name_;
  }


  public int getScans()
  {
    return scans_;
  }


  public float getRtStart()
  {
    return rtStart_;
  }


  public float getRt()
  {
    return rt_;
  }


  public float getRtStop()
  {
    return rtStop_;
  }


  public double getMz()
  {
    return mz_;
  }


  public float getArea()
  {
    return area_;
  }


  public String getAdduct()
  {
    return adduct_;
  }


  public String getIsotope()
  {
    return isotope_;
  }


  public float getScore()
  {
    return score_;
  }


  public float getSignalNoise()
  {
    return signalNoise_;
  }


  public String getDialClassName()
  {
    return dialClassName_;
  }


  public String getDialMs1Name()
  {
    return dialMs1Name_;
  }


  public String getDialMs2Name()
  {
    return dialMs2Name_;
  }


  public String getLdaClassName()
  {
    return ldaClassName_;
  }


  public String getLdaMs1Name()
  {
    return ldaMs1Name_;
  }


  public String getLdaMs2Name()
  {
    return ldaMs2Name_;
  }


  public String getMsmsSpectrum()
  {
    return msmsSpectrum_;
  }

  
  
}
