/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp 
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

package at.tugraz.genome.lda.vos;

import java.util.Hashtable;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class AbsoluteSettingsVO
{
  
  private Hashtable<String,ProbeVolConcVO> volumeSettings_;
  private Hashtable<String, LipidClassSettingVO> classSettings_;
  private Hashtable<String,String> chosenClass_;
  
  public AbsoluteSettingsVO(Hashtable<String,ProbeVolConcVO> volumeSettings_,
      Hashtable<String,LipidClassSettingVO> classSettings_,
      Hashtable<String,String> chosenClass)
  {
    super();
    this.volumeSettings_ = volumeSettings_;
    this.classSettings_ = classSettings_;
    this.chosenClass_ = chosenClass;
  }
  
  public Hashtable<String,ProbeVolConcVO> getVolumeSettings()
  {
    return volumeSettings_;
  }
  public Hashtable<String,LipidClassSettingVO> getClassSettings()
  {
    return classSettings_;
  }
  
  public String getChosenClass(String originalClassName)
  {
  	return chosenClass_.get(originalClassName);
  }
  
}
