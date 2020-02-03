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

import java.util.Vector;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LdaLBLASTCompareVO
{
  private String correctName_;
  private String ldaName_;
  private String lbName_;
  private String rt_;
  private String ldaRts_;
  private double lbProb_;
  private String lbRts_;
  private boolean fp_;
  private int ldaIdentCode_;
  private int lbIdentCode_;
  private boolean ldaMs1Only_;
  private float mz_;
  private boolean fasInOtherCombination_;
  private Vector<LdaLBLASTCompareVO> ms2Evidence_;
  
  public final static int FOUND_NO_FPS = 2;
  public final static int FOUND_BUT_FPS = 1;
  public final static int NOT_FOUND_OR_MS1_ONLY = 0;
  public final static int FALSE_NEGATIVE = -1;
  public final static int FALSE_POSITIVE = -2;
  public final static int FALSE_POSITIVE_AND_NEGATIVE = -3;
  
  private String adduct_;
  private boolean ignoreLDA_;
  private boolean ignoreLB_;

  
  public LdaLBLASTCompareVO(String correctName, String ldaName, String lbName, int ldaIdentCode, int lbIdentCode, String rt, String ldaRts, String lbRts,
      boolean fp, double lbProb, float mz, boolean ldaMs1Only){
    this.rt_ = rt;
    this.correctName_ = correctName;
    this.ldaName_ = ldaName;
    this.lbName_ = lbName;
    this.ldaIdentCode_ = ldaIdentCode;
    this.lbIdentCode_ = lbIdentCode;
    this.fp_ = fp;
    this.ldaRts_ = ldaRts;
    this.lbRts_ = lbRts;
    this.lbProb_ = lbProb;
    this.mz_ = mz;
    this.ldaMs1Only_ = ldaMs1Only;
    ms2Evidence_ = new Vector<LdaLBLASTCompareVO>();
    ignoreLDA_ = false;
    ignoreLB_ = false;
  }

  
  public String getCorrectName()
  {
    return correctName_;
  }

  public boolean isLdaFound()
  {
    return ldaIdentCode_>0;
  }

  public boolean isLbFound()
  {
    return ldaIdentCode_>0;
  }

  public String getLdaName()
  {
    return ldaName_;
  }

  public String getLbName()
  {
    return lbName_;
  }
  

  public int getLdaIdentCode()
  {
    return ldaIdentCode_;
  }


  public int getLbIdentCode()
  {
    return lbIdentCode_;
  }


  public String getRt(){
    return rt_;
  }

  
  public String getLdaRts()
  {
    return ldaRts_;
  }



  public double getLbProb()
  {
    return lbProb_;
  }

  public String getLbRts()
  {
    return lbRts_;
  }
  
public boolean isFp()
  {
    return fp_;
  }


  public void setLbRts(String lbRts)
  {
    this.lbRts_ = lbRts;
  }

  public boolean isLdaMs1Only()
  {
    return ldaMs1Only_;
  }

  public boolean areFasInOtherCombination()
  {
    return fasInOtherCombination_;
  }

  public void setFasInOtherCombination(boolean fasInOtherCombination)
  {
    this.fasInOtherCombination_ = fasInOtherCombination;
  }

  public Vector<LdaLBLASTCompareVO> getMs2Evidence()
  {
    return ms2Evidence_;
  }


  public void setMs2Evidence(Vector<LdaLBLASTCompareVO> ms2Evidence)
  {
    this.ms2Evidence_ = ms2Evidence;
  }


  public boolean isIgnoreLDA()
  {
    return ignoreLDA_;
  }


  public void setIgnoreLDA(boolean ignoreLDA)
  {
    this.ignoreLDA_ = ignoreLDA;
  }


  public boolean isIgnoreLB()
  {
    return ignoreLB_;
  }


  public void setIgnoreLB(boolean ignoreLB)
  {
    this.ignoreLB_ = ignoreLB;
  }


  public String getAdduct()
  {
    return adduct_;
  }


  public void setAdduct(String adduct)
  {
    this.adduct_ = adduct;
  }


  public float getMz()
  {
    return mz_;
  }


  public void setMz(float mz)
  {
    this.mz_ = mz;
  }
  
  
  
  
}
