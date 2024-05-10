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

import java.util.Comparator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class SpectrumPointVO implements Comparable<SpectrumPointVO>
{
  private String mzOriginal_;
  private float mz_;
  private float intensity_;
  
  public SpectrumPointVO(String mzOriginal, float mz)
  {
    super();
    this.mzOriginal_ = mzOriginal;
    this.mz_ = mz;
    this.intensity_ = 0;
  }

  public String getMzOriginal()
  {
    return mzOriginal_;
  }

  public float getMz()
  {
    return mz_;
  }

  public float getIntensity()
  {
    return intensity_;
  }
  
  public void addIntensity(float intensity){
    intensity_ += intensity;
  }

	@Override
	public int compareTo(SpectrumPointVO arg0)
	{
		return Comparator.comparing(SpectrumPointVO::getMz)
				.compare(this, arg0);
	}
  
  
}
