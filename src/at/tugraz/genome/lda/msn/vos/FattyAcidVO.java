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

package at.tugraz.genome.lda.msn.vos;

import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * value object holding information about a lipid acid chain
 * i.e. # C atoms, # double bonds, mass value, chemical formula
 * @author Juergen Hartler
 *
 */
public class FattyAcidVO
{
  // prefix to separate fatty acid
  private String prefix_;
  // number of C atoms
  private int cAtoms_;
  // number of double bonds
  private int doubleBonds_;
  // the mass value
  private double mass_;
  // chemical formula
  private String formula_;
  
  /**
   * all the information has to be provided in the constructor
   * @param cAtoms number of C atoms
   * @param doubleBonds number of double bonds
   * @param mass mass value
   * @param formula chemical formula
   */
  public FattyAcidVO(String prefix, int cAtoms, int doubleBonds, double mass, String formula){
    this.prefix_ = prefix;
    this.cAtoms_ = cAtoms;
    this.mass_ = mass;
    this.doubleBonds_ = doubleBonds;
    this.formula_ = formula;
  }
  
  
  /**
   * 
   * @return the prefix
   */
  public String getPrefix()
  {
    return prefix_;
  }

  /**
   * 
   * @return number of C atoms
   */
  public int getcAtoms()
  {
    return cAtoms_;
  }

  /**
   * 
   * @return number of double bonds
   */
  public int getDoubleBonds()
  {
    return doubleBonds_;
  }

  /**
   * 
   * @return mass value
   */
  public double getMass()
  {
    return mass_;
  }
  
  /**
   * 
   * @return chemical formula
   */
  public String getFormula()
  {
    return formula_;
  }

  /**
   * 
   * @return name consisting of number of C atoms : number of double bonds
   */
  public String getName(){
    return StaticUtils.generateLipidNameString(prefix_+String.valueOf(cAtoms_),doubleBonds_);
  }
  
  /**
   * @return String representing the values of the value object
   */
  public String toString(){
    return ("Name: "+getName()+" Mass: "+mass_+" Formula: "+formula_);
  }
  
}
