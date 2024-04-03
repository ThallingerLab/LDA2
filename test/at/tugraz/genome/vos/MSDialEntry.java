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
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MSDialEntry
{
	
  public final static String MSDIAL_VERSION_4_0 = "4.0";
  public final static String MSDIAL_VERSION_4_9 = "4.9";

  
  private static Hashtable<String,String> adductLookup_ = new Hashtable<String,String>(){{
    put("[M-H]-","-H");
    put("[M+HCOO]-","HCOO");
    put("[M+H]+","H");
    put("[M+Na]+","Na");
    put("[M+H-H2O]+","-OH");
    put("[M+NH4]+","NH4");
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
      float area, String adduct, String isotope, float score, float signalNoise, String msmsSpectrum, String msDialVersion) throws MSDialException
  {
    super();
    this.id_ = id;
    this.name_ = name;
    if (msDialVersion.equalsIgnoreCase(MSDIAL_VERSION_4_0))
    	this.categorizeNameVersion_4_0(mz,rt);
    else if (msDialVersion.equalsIgnoreCase(MSDIAL_VERSION_4_9))
    	this.categorizeNameVersion_4_9(mz,rt);
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
  
  
  private void categorizeNameVersion_4_9(double mz, float rt) throws MSDialException{
    dialClassName_ = null;
    dialMs2Name_ = null;
    ldaClassName_ = null;
    ldaMs1Name_ = null;
    ldaMs2Name_ = null;
    dialMs1Name_ = name_;
    String prefix = null;
    if (dialMs1Name_.startsWith("w/o MS2:"))
      dialMs1Name_ = dialMs1Name_.substring("w/o MS2:".length());
    
    dialClassName_ = dialMs1Name_;
    if (dialClassName_.indexOf(" ")!=-1)
      dialClassName_ = dialClassName_.substring(0,dialClassName_.indexOf(" "));
    if (!dialMs1Name_.contains(" "))
      return;
    if (dialMs1Name_.length()<(dialClassName_.length()+2))
      return;
    char subsequentToEmptySpace = dialMs1Name_.substring(dialClassName_.length()+2,dialClassName_.length()+3).toCharArray()[0];
    //the character subsequent to the empty space must start with a digit for the supported naming convention
    if (!Character.isDigit(subsequentToEmptySpace)) {
      if (dialMs1Name_.length()<=(dialClassName_.length()+4))
      	return;
      prefix = dialMs1Name_.substring(dialClassName_.length()+1,dialClassName_.length()+3);
      if (!prefix.contentEquals("O-")&&!prefix.contentEquals("P-"))
      	return;
      subsequentToEmptySpace = dialMs1Name_.substring(dialClassName_.length()+3,dialClassName_.length()+4).toCharArray()[0];
      if (!Character.isDigit(subsequentToEmptySpace))
      	return;
    }

    //System.out.println(dialMs1Name_);

    //in this case, both the MSDIAL MS1 and MS2 name are stored in the file
    if (dialMs1Name_.contains("|")) {
      dialMs2Name_ = dialMs1Name_.substring(dialMs1Name_.indexOf("|")+1);
      dialMs1Name_ = dialMs1Name_.substring(0,dialMs1Name_.indexOf("|"));
      //System.out.println("1. "+dialMs1Name_+"   "+dialMs2Name_);
    }else {
      //in this case it is assumed that MSDIAL displays an MS2 name only
      if (dialMs1Name_.contains("_")||dialMs1Name_.contains("/")) {
        dialMs2Name_ = dialMs1Name_;
        dialMs1Name_ = contstructMSDIALMs1NameFromMs2Name(dialClassName_, dialMs2Name_, prefix);
        //System.out.println("2. "+dialMs1Name_+"   "+dialMs2Name_);
        //in this case it is assume that MSDIAL MSDIAL displays an MS1 name only
      } else {
        //here, I have to enter if I want to do something specific if no MS2 name is present
        //System.out.println("3. "+dialMs1Name_+"   "+dialMs2Name_);
      }
    }
    dialMs1Name_ = dialMs1Name_.substring(dialClassName_.length()+1);
    if (dialMs2Name_!=null)
      dialMs2Name_ = dialMs2Name_.substring(dialClassName_.length()+1);
    

    
    ldaClassName_ = (prefix!=null ? prefix : "")+dialClassName_;
    ldaMs1Name_ = (prefix!=null ? dialMs1Name_.substring(prefix.length()) : dialMs1Name_);
    if (dialMs1Name_.indexOf(";")>-1) {
      String oh = dialMs1Name_.substring(dialMs1Name_.indexOf(";")+1);
      if (!oh.startsWith("O"))
        throw new MSDialException("For a correct MSDIAL-name, for the first chain, there must be an \"O\" subsequent to the \";\"");
      int nrOh = 1;
      oh = oh.substring(1);
      if (oh.length()>0) {
        try {
          nrOh = Integer.parseInt(oh);
        }catch(NumberFormatException ex) {
          throw new MSDialException("For a correct MSDIAL-name, for the first chain, there must be a number following the \";O\"");
        }
      }
      String hEncoding = "";
      if (nrOh>0) {
        if (nrOh==1)
          hEncoding = "m";
        else if (nrOh==2)
          hEncoding = "d";
        else if (nrOh==3)
          hEncoding = "t";
        else if (nrOh==4)
          hEncoding = "q";
        else if (nrOh==5)
          hEncoding = "p";
        else
          throw new MSDialException("This amount of hydroxylation encodings is not possible in MS-DIAL: "+nrOh+"  "+name_);
      }
      ldaMs1Name_ = hEncoding+dialMs1Name_.substring(0,dialMs1Name_.indexOf(";"));
    }

    

    if (dialMs2Name_!=null){
      String firstFAString = dialMs2Name_;
      String secondFAString = null;
      String thirdFAString = null;
      if (dialMs2Name_.contains("_")||dialMs2Name_.contains("/")) {
        String[] splitted = dialMs2Name_.split("_|/");
        if (splitted.length!=2 && splitted.length!=3)
          throw new MSDialException("This is not a valid MSDIAL-MS2-Name: "+dialMs2Name_);
        firstFAString = splitted[0];
        secondFAString = splitted[1];
        if (splitted.length==3)
        	thirdFAString = splitted[2];
      }
      FattyAcidVO firstFA = getFirstMSDialFattyAcid(firstFAString,dialMs2Name_);
      String hEncoding = "";
      if (firstFA.getOhNumber()>0) {
        if (firstFA.getOhNumber()==1)
          hEncoding = "m";
        else if (firstFA.getOhNumber()==2)
          hEncoding = "d";
        else if (firstFA.getOhNumber()==3)
          hEncoding = "t";
        else if (firstFA.getOhNumber()==4)
          hEncoding = "q";
        else if (firstFA.getOhNumber()==5)
          hEncoding = "p";
        else
        throw new MSDialException("This amount of hydroxylation encodings is not possible in MS-DIAL: "+firstFA.getOhNumber()+"  "+name_);
      }
      
      ldaMs2Name_ = hEncoding+(firstFA.getChainType()==LipidomicsConstants.CHAIN_TYPE_FA_ALKYL ? "O-" : "")+(firstFA.getChainType()==LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL ? "P-" : "")
      		+String.valueOf(firstFA.getcAtoms())+":"+String.valueOf(firstFA.getDoubleBonds());
      if (secondFAString!=null) {
        if (dialMs2Name_.contains("/") || prefix!=null)
          ldaMs2Name_ += "/";
        else
          ldaMs2Name_ += "_";
        FattyAcidVO secondFA = getSecondMSDialFattyAcid(secondFAString,dialMs2Name_);
        if (hEncoding!=null && hEncoding.length()>0)
        	ldaMs2Name_ += (secondFA.getOhNumber()==1 ? "h" : "n");
        ldaMs2Name_ += String.valueOf(secondFA.getcAtoms())+":"+String.valueOf(secondFA.getDoubleBonds());
      }
      if (thirdFAString!=null) {
        if (dialMs2Name_.contains("/"))
          ldaMs2Name_ += "/";
        else
          ldaMs2Name_ += "_";
        FattyAcidVO thirdFA = getSecondMSDialFattyAcid(thirdFAString,dialMs2Name_);
        ldaMs2Name_ += String.valueOf(thirdFA.getcAtoms())+":"+String.valueOf(thirdFA.getDoubleBonds());
      }
    }
    if (ldaClassName_.equalsIgnoreCase("CerP")) {
      ldaClassName_ = "Cer1P";
    } else if (dialClassName_.equalsIgnoreCase("SPB")) {
      ldaClassName_ = "SphBase";
    } else if (ldaClassName_.equalsIgnoreCase("O-PC")) {
      ldaClassName_ = "P-PC";
    }

//    System.out.println(dialClassName_+" ! "+dialMs1Name_+" ! "+dialMs2Name_);
//    System.out.println(ldaClassName_+" ! "+ldaMs1Name_+" ! "+ldaMs2Name_);
    
  }
  
  
  private String contstructMSDIALMs1NameFromMs2Name(String className, String dialMs2Name, String prefix) throws MSDialException {
    String woClass = dialMs2Name.substring(className.length()+1);
    String[] splitted = woClass.split("_|/");
    if (splitted.length!=2 && splitted.length!=3)
      throw new MSDialException("This is not a valid MSDIAL-MS2-Name: "+dialMs2Name);
    //System.out.println(splitted[0]+"         "+splitted[1]);
    int totalOh = 0;
    int totalC = 0;
    int totalDbs = 0;
    String ms1Name = className+" ";
    if (prefix!=null) {
    	ms1Name+=prefix;
    	splitted[0] = splitted[0].substring(prefix.length());
    }
    	
    FattyAcidVO fa1 = getFirstMSDialFattyAcid(splitted[0],dialMs2Name);
    totalC += fa1.getcAtoms();
    totalDbs += fa1.getDoubleBonds();
    totalOh += fa1.getOhNumber();
    
    FattyAcidVO fa2 = getSecondMSDialFattyAcid(splitted[1],dialMs2Name);
    totalC += fa2.getcAtoms();
    totalDbs += fa2.getDoubleBonds();
    totalOh += fa2.getOhNumber();
    
    FattyAcidVO fa3 = null;
    if (splitted.length==3) {
      fa3 = getSecondMSDialFattyAcid(splitted[2],dialMs2Name);
      totalC += fa3.getcAtoms();
      totalDbs += fa3.getDoubleBonds();
      totalOh += fa3.getOhNumber();    	
    }
    
    ms1Name += String.valueOf(totalC)+":"+String.valueOf(totalDbs);
    if (totalOh>0) {
      ms1Name += ";O"+(totalOh>0 ? String.valueOf(totalOh) : "");
    }
    if (fa1.getPrefix()!=null && fa1.getPrefix().startsWith("(d"))
    	ms1Name += fa1.getPrefix();
    if (fa2!=null && fa2.getPrefix()!=null && fa2.getPrefix().startsWith("(d"))
    	ms1Name += fa2.getPrefix();
    if (fa3!=null && fa3.getPrefix()!=null && fa3.getPrefix().startsWith("(d"))
    	ms1Name += fa3.getPrefix();

    return ms1Name;
  }
  
  private FattyAcidVO getFirstMSDialFattyAcid(String firstFA, String fullMS2Name) throws MSDialException {
    int nrOh = 0;
    int nrC = 0;
    int nrDbs = 0;
    short chainType = LipidomicsConstants.CHAIN_TYPE_FA_ACYL;
    if (firstFA.indexOf(";")>-1) {
      String oh = firstFA.substring(firstFA.indexOf(";")+1);
      firstFA = firstFA.substring(0,firstFA.indexOf(";"));
      if (!oh.startsWith("O"))
        throw new MSDialException("For a correct MSDIAL-name, for the first chain, there must be an \"O\" subsequent to the \";\"");
      nrOh = 1;
      oh = oh.substring(1);
      if (oh.length()>0) {
        try {
          nrOh = Integer.parseInt(oh);
        }catch(NumberFormatException ex) {
          throw new MSDialException("For a correct MSDIAL-name, for the first chain, there must be a number following the \";O\"");
        }
      }
    }
    try {
    	String beforeColon = firstFA.substring(0,firstFA.indexOf(":"));
    	if (beforeColon.startsWith("O-")) {
    		beforeColon = beforeColon.substring(2);
    		chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKYL;
    	}
    	if (beforeColon.startsWith("P-")) {
    		beforeColon = beforeColon.substring(2);
    		chainType = LipidomicsConstants.CHAIN_TYPE_FA_ALKENYL;
    	}
      nrC = Integer.parseInt(beforeColon);
      nrDbs = Integer.parseInt(firstFA.substring(firstFA.indexOf(":")+1));
    }catch(NumberFormatException ex) {
      throw new MSDialException("There is something wrong with the first chain in the following MSDIAL-MS2-name: "+fullMS2Name);
    }
    return new FattyAcidVO(chainType, null, nrC, nrDbs, nrOh, -1d, null, null);
  }
  
  private FattyAcidVO getSecondMSDialFattyAcid(String secondFA, String fullMS2Name) throws MSDialException {
    int nrOh = 0;
    int nrC = 0;
    int nrDbs = 0;
    String postColon;
    String prefix = null;

    if (secondFA.indexOf("(2OH)")!=-1) {
      nrOh++;
      secondFA = secondFA.substring(0,secondFA.indexOf("(2OH)"));
    }
    try {
      nrC += Integer.parseInt(secondFA.substring(0,secondFA.indexOf(":")));
      postColon = secondFA.substring(secondFA.indexOf(":")+1);
      if (postColon.contains("(d")) {
      	prefix = postColon.substring(postColon.indexOf("(d"));
      	postColon = postColon.substring(0,postColon.indexOf("(d"));
      }
      	
      nrDbs += Integer.parseInt(postColon);
    }catch(NumberFormatException ex) {
      throw new MSDialException("There is something wrong with the second chain in the following MSDIAL-MS2-name: "+fullMS2Name);
    }
    return new FattyAcidVO((short)-1, prefix, nrC, nrDbs, nrOh, -1d, null, null);
  }
  
  private void categorizeNameVersion_4_0(double mz, float rt) throws MSDialException{
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
