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

package at.tugraz.genome.lda.swing.rdi;

import javax.swing.JFrame;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.PlainDocument;

import at.tugraz.genome.lda.WarningMessage;

/**
 * is for a JTextField for entering integer values that are between a specifiable range
 * @author Juergen Hartler
 *
 */
public class IntegerRangeDocument extends PlainDocument
{

  private static final long serialVersionUID = -6496137843055033710L;
  
  /** the lowest allowed integer value*/
  private int min_;
  /** the highest allowed integer value*/
  private int max_;
  
  /**
   * constructor for setting the specifiable integer range
   * @param min lowest allowed integer value
   * @param max highest allowed integer value
   */
  public IntegerRangeDocument(int min, int max){
    this.min_ = min;
    this.max_ = max;
  }
  
  public void insertString(int offs, String str, AttributeSet a) throws BadLocationException{
    String newText = getCleanedText(str);
    if (newText.equalsIgnoreCase("otherChars")) return;
    String oldText = this.getText(0, this.getLength());
    String textResult = oldText.substring(0,offs)+newText+oldText.substring(offs);
    if (checkInput(textResult,false)) super.insertString(offs, newText, a);
  }
  
  public void remove(int offs, int length) throws BadLocationException{
    String oldText = this.getText(0, this.getLength());
    String textResult =  oldText.substring(0,offs)+oldText.substring(offs+length);
    if (checkInput(textResult,true)) super.remove(offs, length);
  }
  
  public void replace(int offs, int length, String str, AttributeSet a) throws BadLocationException{
    String newText = getCleanedText(str);
    if (newText.equalsIgnoreCase("otherChars")) return;
    String oldText = this.getText(0, this.getLength());
    String textResult = oldText.substring(0,offs)+newText+oldText.substring(offs+length);
    if (checkInput(textResult,false)) super.insertString(offs, newText, a);

  }
  
  /**
   * cleans any non-integer characters from the input string
   * @param str the input text
   * @return the cleaned String
   */
  private String getCleanedText(String str){
    char[] chars = str.trim().toCharArray();
    String newText = "";
    boolean otherChars = false;
    for (int i=0; i!=chars.length; i++){
      if (Character.isDigit(chars[i])) newText += chars[i];
      else otherChars = true;
    }
    if (otherChars){
      new WarningMessage(new JFrame(), "Error", "Only digits are allowed for this field");
      return "otherChars";
    }
    return newText;
  }
  
  /**
   * checks if the final inserted value is integer and is inside the allowed range
   * @param textResult the text that will be in the field after the insert/remove/replace operation
   * @param allowNull is a null value or empty string allowed
   * @return true if the input value is integer and is inside the allowed range
   */
  private boolean checkInput(String textResult, boolean allowNull){
    if (textResult==null || textResult.length()==0) return allowNull;
    try{
      int number = Integer.parseInt(textResult);
      if (number<min_){
        new WarningMessage(new JFrame(), "Error", "The value must not be smaller than "+min_);
        return false;
      }
      if (number>max_){
        new WarningMessage(new JFrame(), "Error", "The value must not be higher than "+max_);
        return false;
      }
      return true;
    } catch (NumberFormatException nfx){
      return false;
    }
  }
}
