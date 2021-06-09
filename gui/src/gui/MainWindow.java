package gui;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.border.BevelBorder;
import javax.swing.border.LineBorder;

import lcip.LCIP;
//import lcip.Rip;

import java.awt.Color;

/*
 * Authors: José Duarte 	(hfduarte@ua.pt)
 * 			Mark Mckenney 	(marmcke@siue.edu)
 * 
 * Version: 1.0 (2020-05-29)
 * 
 */

public class MainWindow extends JFrame 
{
	private static final long serialVersionUID = -7003878987581378628L;

	private JPanel contentPane;
	private JFrame f;
	private JButton lcip_btn;
	private JButton rip_btn;
	
	public MainWindow() 
	{
		f = this;
		
		draw_UI();
		add_listeners();
		
		// load configuration.
		
		//   ->
	}
	
	public void draw_UI() 
	{
		setTitle("Lines and Concavities Interpolation Problem");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    setSize(830, 700);
	    setLocationRelativeTo(null);
	    setExtendedState(f.getExtendedState() | JFrame.MAXIMIZED_BOTH);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);
		contentPane.setLayout(new GridLayout(8, 1, 10, 10));
				
		JPanel format = new JPanel();
		format.setBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		contentPane.add(format);

		JPanel at_p = new JPanel();
		at_p.setBorder(new LineBorder(new Color(0, 0, 0)));
		contentPane.add(at_p);
		
		lcip_btn = new JButton("LCIP");
		at_p.add(lcip_btn);
		
		rip_btn = new JButton("RIP");
		rip_btn.setEnabled(false);
		at_p.add(rip_btn);
	}

	public void add_listeners()
	{		
		lcip_btn.addActionListener(new ActionListener() 
		{
            public void actionPerformed(ActionEvent event) 
            {   
        		LCIP lcip = new LCIP();
        		lcip.setVisible(true);
            }
        });
		
		/*
		rip_btn.addActionListener(new ActionListener() 
		{
            public void actionPerformed(ActionEvent event) 
            {   
        		Rip rip = new Rip();
        		rip.setVisible(true);
            }
        });
        */
	}
}