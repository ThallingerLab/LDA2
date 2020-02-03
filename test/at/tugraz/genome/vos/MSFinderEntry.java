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

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MSFinderEntry
{
  
  private double rt_ = -1d;
  private String title_ = null;
  private int ms1Count_ = -1;
  private int msmsCount_ = -1;
  private float mz_ = -1f;
  private String adduct_ = null;
  
  private Hashtable<Integer,MSFinderHitVO> hits_;
  
    
  public MSFinderEntry(String filePath, String title, int ms1Count, int msmsCount, float mz, String adduct)
  {
    super();
    String path = filePath;
    int indSlash = path.lastIndexOf("/");
    int indBSlash = path.lastIndexOf("\\");
    if (indBSlash>indSlash)
      indSlash = indBSlash;
    path = path.substring(indSlash+1);
    path = path.substring(path.indexOf("_")+1);
    rt_ = Double.parseDouble(path.substring(0,path.indexOf("_")));
    this.title_ = title;
    this.ms1Count_ = ms1Count;
    this.msmsCount_ = msmsCount;
    this.mz_ = mz;
    this.adduct_ = adduct;
    this.hits_ = new Hashtable<Integer,MSFinderHitVO>();
  }
  
  
  public double getRt()
  {
    return rt_;
  }
  public String getTitle()
  {
    return title_;
  }
  public int getMs1Count()
  {
    return ms1Count_;
  }
  public int getMsmsCount()
  {
    return msmsCount_;
  }
  public float getMz()
  {
    return mz_;
  }
  public String getAdduct()
  {
    return adduct_;
  }
  public void addHit(MSFinderHitVO hitVO) {
    this.hits_.put(hitVO.getRank(), hitVO);
  }


  public Hashtable<Integer,MSFinderHitVO> getHits()
  {
    return hits_;
  }
  
}
