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

package at.tugraz.genome.lda.swing;

import java.awt.Color;
import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JButton;
import javax.swing.border.LineBorder;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ExportButton extends JButton implements MouseListener
{
  private static final long serialVersionUID = 211052889655909533L;

  public ExportButton(String name,String actionCommand, Color fontColor, Color bkgrnd, ActionListener parent) {
    super(name);
    setForeground(fontColor);
    setActionCommand(actionCommand);
    addActionListener(parent);
    Font textFont = new Font("Helvetica",Font.BOLD, 9);
    //setFont(getFont().deriveFont(10f));
    setFont(textFont);
    setMargin(new Insets(1,5,1,5));
    setBackground(bkgrnd);
    setMargin(new Insets(1,15,1,15));

    this.setBorder(new LineBorder(fontColor));
    setBorderPainted(false);
    this.setRolloverEnabled(false);
    this.setContentAreaFilled(false);
    this.setOpaque(false);
    this.addMouseListener(this);
    this.setFocusPainted(false);
    this.setContentAreaFilled(false);
  }
  
  public void mouseEntered(MouseEvent e) {
    this.setText("<html><u>"+getText()+"</u></html>");
    //setBorderPainted(true);
  }


  public void mouseExited(MouseEvent e) {
    this.setText(getText().substring("<html><u>".length(),getText().length()-"</u></html>".length()));
  }

  public void mouseClicked(MouseEvent e)
  {
  }

  public void mousePressed(MouseEvent e)
  {
  }

  public void mouseReleased(MouseEvent e)
  {
  }

  
}
