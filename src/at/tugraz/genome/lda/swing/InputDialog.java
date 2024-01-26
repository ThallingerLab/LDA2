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

package at.tugraz.genome.lda.swing;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import at.tugraz.genome.lda.TooltipTexts;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class InputDialog extends JDialog implements ActionListener
{
  
  private static final long serialVersionUID = 8515575443848106452L;
  
  private String originalValue_;
  JTextField text_;
  private String enteredValue_;

  public InputDialog(JFrame parent, String title, String message, String value) {
    super(parent, title, true);
/*    if (parent != null) {
      System.out.println("HAllo");
      Dimension parentSize = parent.getSize(); 
      Point p = parent.getLocation(); 
      setLocation(p.x + parentSize.width / 2, p.y + parentSize.height / 2);
    }*/
    originalValue_ = value;
    enteredValue_ = value;
    setLocation(380,240);
    JPanel messagePane = new JPanel();
    messagePane.add(new JLabel(message));
    getContentPane().add(messagePane,BorderLayout.NORTH);
    JPanel inputPane = new JPanel();
    text_ = new JTextField(15);
    text_.setText(value);
    inputPane.add(text_); 
    getContentPane().add(inputPane, BorderLayout.CENTER);
    JButton button = new JButton("OK"); 
    JPanel buttonPane = new JPanel();
    button.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    buttonPane.add(button); 
    button.addActionListener(this);
    getContentPane().add(buttonPane, BorderLayout.SOUTH);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    pack(); 
    setVisible(true);
  }  
  
  public void actionPerformed(ActionEvent e) {
    this.enteredValue_ = text_.getText();
    setVisible(false); 
    dispose(); 
  }
  
  public String getEnteredText(){
    if (enteredValue_!=null && enteredValue_.length()>0){
      return enteredValue_;
    }else
      return originalValue_;
  }
}
