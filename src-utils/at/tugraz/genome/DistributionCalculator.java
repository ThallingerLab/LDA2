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
package at.tugraz.genome;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.parser.exceptions.SpectrummillParserException;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class DistributionCalculator extends JFrame implements ActionListener
{


  private static final long serialVersionUID = -336459981982138449L;

  private JPanel displayPanel_ = new JPanel();
  private JTextField input_;

  public DistributionCalculator(){
    JFrame frame = new JFrame("Distribution Calculator");
    frame.setSize(600,300);
    frame.setLocationByPlatform(true);
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    JPanel mainWindow = new JPanel();
    mainWindow.setLayout(new BorderLayout());
    
    JPanel inputPanel = new JPanel();
    JLabel label = new JLabel("Enter chem. formula: ");
    inputPanel.add(label);
    input_ = new JTextField(20);
    inputPanel.add(input_);
    JButton button = new JButton("OK");
    button.addActionListener(this);
    button.setActionCommand("calc");
    button.setFont(button.getFont().deriveFont(10f));
    button.setMargin(new Insets(1,5,1,5));
    inputPanel.add(button);
    
    displayPanel_ = new JPanel();
    displayPanel_.setLayout(new GridBagLayout());
    mainWindow.add(displayPanel_,BorderLayout.CENTER);
    mainWindow.add(inputPanel,BorderLayout.NORTH);
    
    frame.setVisible(true);
    frame.add(mainWindow);
  }
  
  
  /**
   * @param args
   */
  public static void main(String[] args)
  {
    new DistributionCalculator();
  }


  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("calc")){
      displayDistribution();
    }
    
  }

  private void displayDistribution(){
    if (input_.getText()==null || input_.getText().length()==0){
      new WarningMessage(new JFrame(), "Error", "Please enter a formula!");
      return;
    }
    String formula = input_.getText();
    Vector<Vector<Double>> bothDistris = calculateDistribution(formula);
    displayPanel_.removeAll();
    Vector<Double> distri = bothDistris.get(0);
    for (int i=0; i!=distri.size(); i++){
      JLabel label = new JLabel("+"+i+": ");
      JLabel value = new JLabel(distri.get(i).toString());
      displayPanel_.add(label,new GridBagConstraints(0, i, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 2, 0));
      displayPanel_.add(value,new GridBagConstraints(1, i, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 2, 0));
    }
    int add = distri.size()+1;
    Vector<Double> mass = bothDistris.get(1);
    if (bothDistris.size()>2){
      distri = bothDistris.get(1);
      mass = bothDistris.get(2);
      for (int i=0; i!=distri.size(); i++){
        JLabel label = new JLabel("-"+i+": ");
        JLabel value = new JLabel(distri.get(i).toString());
        displayPanel_.add(label,new GridBagConstraints(0, i+add, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 2, 0));
        displayPanel_.add(value,new GridBagConstraints(1, i+add, 1, 1, 0.0, 0.0
          ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 2, 0));
      }
      add += distri.size();
    }
    displayPanel_.add(new JLabel("mass: "),new GridBagConstraints(0, add+1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 2, 0));
      displayPanel_.add(new JLabel(Calculator.FormatNumberToString(mass.get(0),5d)),new GridBagConstraints(1, add+1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 2, 0));
    
    displayPanel_.invalidate();
    displayPanel_.updateUI();
  }
  
  private  Vector<Vector<Double>> calculateDistribution(String formula) {
    Vector<Vector<Double>> distri = new Vector<Vector<Double>>();
    ElementConfigParser aaParser = new ElementConfigParser("elementconfig.xml");
    try {
      Hashtable<String,Integer> categorized = StaticUtils.categorizeFormula(formula);
      String formulaToCheck = "";
      for (String elem : categorized.keySet()){
        if (formulaToCheck.length()>0) formulaToCheck += " ";
        formulaToCheck += elem+categorized.get(elem).toString();
      }
      aaParser.parse();
      distri = aaParser.calculateChemicalFormulaIntensityDistribution(formulaToCheck,5,false);
      Vector<Double> mass = new Vector<Double>();
      mass.add(aaParser.calculateTheoreticalMass(formulaToCheck, false));
      distri.add(mass);
      return distri;
    }
    catch (SpectrummillParserException e) {
      new WarningMessage(new JFrame(), "Warning", e.getMessage());
    }
    catch (ChemicalFormulaException e) {
      new WarningMessage(new JFrame(), "Warning", e.getMessage());
    }
    return distri;
  }
}
