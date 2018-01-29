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

package at.tugraz.genome.lda;

import at.tugraz.genome.lda.xml.RawToChromTranslator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MzxmlToChromThread extends Thread
{
  private String filePath_;
  private boolean finished_ = false;
  private String errorString_;
  private int numberOfThreads_;
  /** true when the mzXML file contains polarity switched data*/
  private boolean polaritySwitched_;
  
  public MzxmlToChromThread(String filePath, int numberOfThreads){
    this.filePath_ = filePath;
    this.numberOfThreads_ = numberOfThreads;
  }
  
  public void run(){
    try{
      polaritySwitched_ = false;
      translateToChrom(filePath_,numberOfThreads_);
    } catch (Exception ex){
      ex.printStackTrace();
      errorString_ = ex.toString();     
    }
    finished_ = true;
  }
  
  public boolean finished(){
    return this.finished_;
  }
  
  public String getErrorString(){
    return this.errorString_;
  }
  
  private void translateToChrom(String filePath, int numberOfThreads) throws Exception{
    RawToChromTranslator translator = new RawToChromTranslator(filePath,"mzXML", LipidomicsConstants.getmMaxFileSizeForChromTranslationAtOnceInMB(),
        numberOfThreads,LipidomicsConstants.getChromMultiplicationFactorForInt(),LipidomicsConstants.getChromLowestResolution(),LipidomicsConstants.isMS2());
      translator.translateToChromatograms();
      polaritySwitched_ = translator.isPolaritySwitched();
  }

  /**
   * 
   * @return true when the mzXML file contains polarity switched data
   */
  public boolean isPolaritySwitched()
  {
    return polaritySwitched_;
  }
  
  
  
}
