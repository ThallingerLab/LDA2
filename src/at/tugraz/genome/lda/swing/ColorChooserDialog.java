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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JPanel;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.interfaces.ColorChangeListener;
import at.tugraz.genome.lda.utils.ColorSequence;
import at.tugraz.genome.maspectras.utils.ColorEncoder;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ColorChooserDialog extends JDialog implements ActionListener
{

  private static final long serialVersionUID = -5108488715688689425L;
  
  public static final int DEFAULT_TYPE = 0;
  public static final int EXPERIMENT_TYPE = 1;
  public static final int GROUP_TYPE = 2;
//  public static final int CLASS_TYPE = 3;
  
  private JColorChooser colorChooser_;
  
  private Color defaultColor_;
  private Hashtable<String,Color> expColors_;
  private Hashtable<String,Color> groupColors_;
//  private Hashtable<String,Color> classesColors_;
  private ColorChangeListener changeListener_;
//  private Hashtable<String,String> fromShortToExpName_;
  
//  SampleLookup sampleLookup_;
  
  private int currentType_;
  private String currentName_;

  public ColorChooserDialog(JFrame parent, String title, Vector<String> expNames, Vector<String> groupNames,
      ColorChangeListener changeListener){
    super(parent, title, true);
    
//    sampleLookup_ = sampleLookup;
    defaultColor_ = Color.WHITE;
    expColors_ = createColorScheme(expNames);
//    this.getDisplayNames();
    groupColors_ = createColorScheme(groupNames);
    changeListener_ = changeListener;
    
    setLocation(380,240);
    this.setLayout(new BorderLayout());
    
    colorChooser_ = new JColorChooser();
    add(colorChooser_,BorderLayout.CENTER);
    
    JPanel buttonPane = new JPanel();
    JButton cancelButton = new JButton("Cancel");
    buttonPane.add(cancelButton); 
    cancelButton.addActionListener(this);
    cancelButton.setActionCommand("CancelColorSelection");
    cancelButton.setToolTipText(TooltipTexts.CANCEL_GENERAL);
    JButton okButton = new JButton("OK"); 
    buttonPane.add(okButton); 
    okButton.addActionListener(this);
    okButton.setActionCommand("AcceptColorSelection");
    okButton.setToolTipText(TooltipTexts.ACCEPT_GENERAL);
    add(buttonPane, BorderLayout.SOUTH);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    pack(); 
    setVisible(false);
  }
  
//  private Vector<String> getDisplayNames(){
//    Vector<String> names = new Vector<String>();
//    fromShortToExpName_ = new Hashtable<String,String>();
//    for (String name : expColors_.keySet()){
//      String displayName = sampleLookup_.getDisplayName(name);
//      names.add(displayName);
//      fromShortToExpName_.put(displayName, name);
//    }
//    return names;
//  }
  
  public void showColorChooser(String title, int type, String name){
    currentType_ = type;
    currentName_ = name;
    setTitle(title);
    Color initColor = getColor(type,name);
    colorChooser_.setColor(initColor);
    setVisible(true);
  }
  
  public void actionPerformed(ActionEvent e)
  {
    if (e.getActionCommand().equalsIgnoreCase("AcceptColorSelection")){
      Color selectedColor = colorChooser_.getColor();
      if (currentType_==DEFAULT_TYPE)
        defaultColor_ = selectedColor;
      else if (currentType_==EXPERIMENT_TYPE)
        expColors_.put(currentName_, selectedColor);
      else if (currentType_==GROUP_TYPE)
        groupColors_.put(currentName_, selectedColor);
      setVisible(false);
      changeListener_.updateColorChange();
    }else if(e.getActionCommand().equalsIgnoreCase("CancelColorSelection")){
      setVisible(false);
    }
  }
  
  @SuppressWarnings("unchecked")
  private Hashtable<String,Color> createColorScheme(Vector<String> valueNames){
    Hashtable<String,Color> colors = new Hashtable<String,Color>();
    if (valueNames!=null){
      Vector<int[]> colorChannels = (new ColorEncoder(valueNames.size())).getChannels();
      Color[] defaultColors = ColorSequence.getDefaultColors();
    
      int count = 0;
      for (String name: valueNames){
        Color color = null;
        if (valueNames.size()>defaultColors.length){
          int[] rgb = colorChannels.get(count);
          color = new Color(rgb[0],rgb[1],rgb[2]);
        }else{
          color = defaultColors[count];
        }
        colors.put(name, color);
        count++;
      }
    }
    return colors;
  }
  
  public Color getColor(int type, String name){
    if (type==DEFAULT_TYPE)
      return defaultColor_;
    else if (type==EXPERIMENT_TYPE)
      return expColors_.get(name);
    else if (type==GROUP_TYPE)
      return groupColors_.get(name);
    else return Color.BLACK;
  }
  
  
}
