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

package at.tugraz.genome.lda.exception;


/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class NumberOutOfRangeException extends Exception
{
  private static final long serialVersionUID = 1L;


  /**
   *  Constructor for NumberOutOfRangeException.
   */
  public NumberOutOfRangeException() {
      super();
  }


  /**
   *  Constructor for NumberOutOfRangeException.
   *
   *@param  message
   */
  public NumberOutOfRangeException(String message) {
      super(message);
  }


  /**
   *  Constructor for NumberOutOfRangeException.
   *
   *@param  message
   *@param  cause
   */
  public NumberOutOfRangeException(String message, Throwable cause) {
      super(message, cause);
  }


  /**
   *  Constructor for NumberOutOfRangeException.
   *
   *@param  cause
   */
  public NumberOutOfRangeException(Throwable cause) {
      super(cause);
  }

}
