/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2024 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

import java.util.Collections;
import java.util.Vector;

import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.utils.StaticUtils;

public class CompareTest
{
	public static void main(String[] args)
  {
		try
		{
//			Vector<FattyAcidVO> vos = StaticUtils.decodeFAsFromHumanReadableName("18:1(n-7)_O-18:1", new HydroxyEncoding("hydroxylationEncoding.txt"),
//					new HydroxyEncoding("lcbHydroxylationEncoding.txt"), false, null);
			Vector<FattyAcidVO> vos = StaticUtils.decodeFAsFromHumanReadableName("18:0_O-18:1(n-7)", new HydroxyEncoding("hydroxylationEncoding.txt"),
					new HydroxyEncoding("lcbHydroxylationEncoding.txt"), false, null);
			
			System.out.println("before");
			for (FattyAcidVO vo : vos)
			{
				System.out.println(vo.toString());
			}
			Collections.sort(vos);
			System.out.println("after");
			for (FattyAcidVO vo : vos)
			{
				System.out.println(vo.toString());
			}
		}
		catch (Exception ex)
		{
			System.out.println("hi");
		}
  }
	
	
}
