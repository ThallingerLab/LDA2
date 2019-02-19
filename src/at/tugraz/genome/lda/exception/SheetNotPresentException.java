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

package at.tugraz.genome.lda.exception;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class SheetNotPresentException extends Exception
{
  
  /**
   * 
   */
  private static final long serialVersionUID = -5180022768477876977L;


  /**
   *  Constructor for SheetNotPresentException.
   */
  public SheetNotPresentException() {
      super();
  }


  /**
   *  Constructor for SheetNotPresentException.
   *
   *@param  message
   */
  public SheetNotPresentException(String message) {
      super(message);
  }


  /**
   *  Constructor for SheetNotPresentException.
   *
   *@param  message
   *@param  cause
   */
  public SheetNotPresentException(String message, Throwable cause) {
      super(message, cause);
  }


  /**
   *  Constructor for SheetNotPresentException.
   *
   *@param  cause
   */
  public SheetNotPresentException(Throwable cause) {
      super(cause);
  }

}
