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

/**
 * 
 * @author Juergen Hartler
 *
 */
public class MSFinderHitVO
{

  private int rank_;
  private String structureName_;
  private float score_;
  private String database_;
  private String formula_;
  private String ontology_;
  private String inchikey_;
  private String smile_;
  
  
  public MSFinderHitVO(int rank, String structureName, float score, String database,
      String formula, String ontology, String inchikey, String smile){
    super();
    this.rank_ = rank;
    this.structureName_ = structureName;
    this.score_ = score;
    this.database_ = database;
    this.formula_ = formula;
    this.ontology_ = ontology;
    this.inchikey_ = inchikey;
    this.smile_ = smile;
  }


  public int getRank()
  {
    return rank_;
  }


  public String getStructureName()
  {
    return structureName_;
  }


  public float getScore()
  {
    return score_;
  }


  public String getDatabase()
  {
    return database_;
  }


  public String getFormula()
  {
    return formula_;
  }


  public String getOntology()
  {
    return ontology_;
  }


  public String getInchikey()
  {
    return inchikey_;
  }


  public String getSmile()
  {
    return smile_;
  }
  
  
  
  
}
