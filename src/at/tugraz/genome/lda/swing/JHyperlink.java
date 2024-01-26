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

import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Desktop;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;


import javax.swing.JLabel;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class JHyperlink extends JLabel implements MouseListener
{

  /**
   * 
   */
  private static final long serialVersionUID = -7224370255158950845L;
  private URI uri = null;
  private Component parentComponent = this;

  /**
   * Constructor
   *
   * @param String text
   * @param String uri
   */

  public JHyperlink( String text, String uri ) {
  this.setText( "<html><u>" + text + "</u></html>" );

  try {

  this.uri = new URI( uri );

  } catch (URISyntaxException e) {

  e.printStackTrace();

  }

  initialize();

  } // end constructor

  /**
   * This method returns the URI
   *
   * @return uri
   */

  public URI getURI() {

  return uri;

  } // end method getURI

  /**
   * This method initializes the
   * link.
   */

  private void initialize() {

  this.setForeground( Color.BLUE );
  this.setOpaque( false );
  this.setToolTipText( getURI().toString() );
  this.addMouseListener( this );

  } // end method initialize

  public void mouseClicked(MouseEvent e) {

  if ( e.getClickCount() > 0 ) {
  if ( Desktop.isDesktopSupported() ) {
  Desktop desktop = Desktop.getDesktop();

  try {
  desktop.browse( getURI() );
  } catch ( IOException ex ) {
  ex.printStackTrace();
  }
  }
  }

  } // end method mouseClicked

  public void mouseEntered(MouseEvent arg0) {
  parentComponent.setCursor( Cursor.getPredefinedCursor( Cursor.HAND_CURSOR ) );
  } // end method mouseEntered


  public void mouseExited(MouseEvent arg0) {
  parentComponent.setCursor( Cursor.getPredefinedCursor( Cursor.DEFAULT_CURSOR ) );
  } // end method mouseExited

  public void mousePressed(MouseEvent arg0) {

  this.setForeground( Color.RED );
  parentComponent.setCursor( Cursor.getPredefinedCursor( Cursor.WAIT_CURSOR ) );

  } // end method mousePressed

  public void mouseReleased(MouseEvent arg0) {

  this.setForeground( Color.BLUE );
  parentComponent.setCursor( Cursor.getPredefinedCursor( Cursor.DEFAULT_CURSOR ) );

  } // end method mouseReleased

  public void setParentComponent(Component parent){
  this.parentComponent = parent;
  }
} // end class JHyperlink


