package tools;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.text.DecimalFormat;
import java.util.ArrayList;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import javax.swing.JTable;
import javax.swing.SwingConstants;

import java.awt.BorderLayout;
import javax.swing.border.LineBorder;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableModel;
import javax.swing.JButton;
import javax.swing.JTextField;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.awt.FlowLayout;

/*
 * Authors: José Duarte 	(hfduarte@ua.pt)
 * 			Mark Mckenney 	(marmcke@siue.edu)
 * 
 * Version: 1.0 (2020-05-29)
 * 
 */

public class Table extends JFrame
{
	private static final long serialVersionUID = 1L;
	
	private String title;
	private JTable table;
	private MyTableModel model;
	
	private String[] col_names;
	private ArrayList<String> data;
	private JTextField csv_path_tx;
	
	private int height;
	private int weight;
	
	public Table(
			String title,
			String[] col_names,
			ArrayList<String> data,
			int weight,
			int height
		) 
	{
		this.title = title;
		this.col_names = col_names;
		this.data = data;
			
		table = new JTable();
		table.setBorder(new LineBorder(new Color(0, 0, 0)));
		getContentPane().add(table, BorderLayout.CENTER);
		
		this.weight = weight;
		this.height = height;

		initUI();
    }

    private void initUI() 
    {   
		setTitle(title);
	    setSize(1024, 720);
	    setLocationRelativeTo(null);
	    setExtendedState(getExtendedState() | JFrame.MAXIMIZED_BOTH);
    	
    	JPanel content = new JPanel();
    	content.setLayout(new BorderLayout(0, 0));
        model = new MyTableModel(data.size(), col_names);
        
        int col_id = 0;
        String row_data = null;
        String[] values = null;
        
        for (int row = 0; row < data.size(); row++) 
        {
        	row_data = data.get(row);
        	values = row_data.split(";");
        	
        	col_id = 0;
        	
        	for(int k = 0; k < values.length; k++)
        	{
        		if(k == values.length - 1)
        			model.setValueAt(values[k], row, col_id);
        		else
        			model.setValueAt(Double.parseDouble(values[k]), row, col_id);
        		
        		/*
        		if(k < 1)
        			model.setValueAt(values[k], row, col_id);
        		else
        			model.setValueAt(Double.parseDouble(values[k]), row, col_id);
        		*/
        		col_id++;
        	}
        }
                
        TableCellRenderer number_renderer = new MyNumberCellRenderer(4);
        TableCellRenderer text_renderer = new MyTextCellRenderer();
        
        JTable table_1 = new JTable(model)
		{
			private static final long serialVersionUID = 1L;

			public Component prepareRenderer(TableCellRenderer renderer, int row, int column)
			{
				Component c = super.prepareRenderer(renderer, row, column);

				//  Alternate row color

				if (!isRowSelected(row))
					c.setBackground(row % 2 == 0 ? new Color(240, 255, 240) : Color.white);
				
				return c;
			}
		};
        
        for(int h = 0; h < col_names.length - 1; h++)
        	table_1.getColumnModel().getColumn(h).setCellRenderer(number_renderer);
        
        table_1.getColumnModel().getColumn(col_names.length - 1).setCellRenderer(text_renderer);
        
        JPanel panel = new JPanel();
        content.add(panel, BorderLayout.SOUTH);
        panel.setLayout(new FlowLayout(FlowLayout.CENTER, 25, 25));
        
        JPanel panel_c = new JPanel();
        content.add(panel_c, BorderLayout.CENTER);
        panel_c.setLayout(new FlowLayout(FlowLayout.CENTER, 25, 25));
        
        JScrollPane scrollPane = new JScrollPane(table_1);
        scrollPane.setPreferredSize(new Dimension(weight, height));
        panel_c.add(scrollPane);
        getContentPane().add(content);
                
        csv_path_tx = new JTextField();
        csv_path_tx.setBackground(new Color(240, 255, 240));
        csv_path_tx.setText("/home/user/output.csv");
        panel.add(csv_path_tx);
        csv_path_tx.setColumns(20);
        
        JButton save_to_csv_bt = new JButton("Save as CSV");
        save_to_csv_bt.setFont(new Font("Tahoma", Font.BOLD, 11));
        save_to_csv_bt.setForeground(new Color(25, 25, 112));
        panel.add(save_to_csv_bt);

        save_to_csv_bt.addActionListener(new ActionListener() 
		{
			public void actionPerformed(ActionEvent arg0) 
		    {
				String filename = csv_path_tx.getText();
				PrintWriter writer = null;
				String column_names = "";
				
				try {
					writer = new PrintWriter(filename, "UTF-8");
				} catch (FileNotFoundException | UnsupportedEncodingException e) 
				{
					JOptionPane.showMessageDialog(null, e.getMessage(), "Save", JOptionPane.INFORMATION_MESSAGE);
					return;
				}
				
				//.
				
				for(String s : col_names)
					column_names += s + ";";
				
				writer.println(column_names);
				
				//.
				
				for(String values : data)
					writer.println(values);
				
				//.
				
				writer.close();
				JOptionPane.showMessageDialog(null, "Saved to " + filename + "!", "Save", JOptionPane.INFORMATION_MESSAGE);	
		    }
		});
    }
    
    public class MyNumberCellRenderer extends DefaultTableCellRenderer 
    {
		private static final long serialVersionUID = 1L;

		private final DecimalFormat format;
		
		/**
         * Default constructor - builds a renderer that right justifies the 
         * contents of a table cell.
         */
        public MyNumberCellRenderer(int decimalPlaces) 
        {
            super();
            setHorizontalAlignment(SwingConstants.RIGHT);
            
    		setOpaque(true);

    		StringBuilder formatString = new StringBuilder();
    		formatString.append("###############0");

    		if (decimalPlaces > 0) {
    			formatString.append(".");
    		}

    		for (int i = 0; i < decimalPlaces; i++) {
    			formatString.append("0");
    		}

    		format = new DecimalFormat(formatString.toString());
        }
        
        /**
         * Returns itself as the renderer. Supports the TableCellRenderer interface.
         *
         * @param table  the table.
         * @param value  the data to be rendered.
         * @param isSelected  a boolean that indicates whether or not the cell is 
         *                    selected.
         * @param hasFocus  a boolean that indicates whether or not the cell has 
         *                  the focus.
         * @param row  the (zero-based) row index.
         * @param column  the (zero-based) column index.
         *
         * @return the component that can render the contents of the cell.
         */
        public Component getTableCellRendererComponent(final JTable table, 
                final Object value, final boolean isSelected, 
                final boolean hasFocus, final int row, final int column) 
        {
    		if (isSelected) {
    			setBackground(table.getSelectionBackground());
    		}
            else 
            {
                //setBackground(null);
            }

            if (value != null) {
                //setText(format.format(value));
                try
                {
                	setText(format.format(value));
                }
                catch(IllegalArgumentException e)
                {
                	setText((String) value);
                }
            }
            else 
            {
                setText("");
            }

            /*
            try
            {
            	setText(format.format(value));
            }
            catch(IllegalArgumentException e)
            {
            	setText((String) value);
            }
            */
    		//setText(format.format(value));
            return this;
        }
    }
    
    public class MyTextCellRenderer extends DefaultTableCellRenderer 
    {
		private static final long serialVersionUID = 1L;
		
		/**
         * Default constructor - builds a renderer that right justifies the 
         * contents of a table cell.
         */
        public MyTextCellRenderer() 
        {
            super();
            setHorizontalAlignment(SwingConstants.CENTER);
            
    		setOpaque(true);
        }
        
        /**
         * Returns itself as the renderer. Supports the TableCellRenderer interface.
         *
         * @param table  the table.
         * @param value  the data to be rendered.
         * @param isSelected  a boolean that indicates whether or not the cell is 
         *                    selected.
         * @param hasFocus  a boolean that indicates whether or not the cell has 
         *                  the focus.
         * @param row  the (zero-based) row index.
         * @param column  the (zero-based) column index.
         *
         * @return the component that can render the contents of the cell.
         */
        public Component getTableCellRendererComponent(final JTable table, 
                final Object value, final boolean isSelected, 
                final boolean hasFocus, final int row, final int column) 
        {
    		if (isSelected) {
    			setBackground(table.getSelectionBackground());
    		}
            else 
            {
                //setBackground(null);
            }
    		
            if (value != null) 
            {
            	setText((String) value);
            }
            else 
            {
                setText("");
            }

            return this;
        }
    }
    
    /**
     * A table model.
     */
    static class MyTableModel extends AbstractTableModel implements TableModel 
    {
		private static final long serialVersionUID = 1L;
		private Object[][] data;
		private String[] col_names;

        /**
         * Creates a new table model
         *
         * @param rows  the row count.
         */
        public MyTableModel(int rows, String[] col_names) 
        {
        	this.col_names = col_names;
            this.data = new Object[rows][col_names.length];
        }

        /**
         * Returns the column count.
         *
         * @return 7.
         */
        public int getColumnCount() {
            return col_names.length;
        }

        /**
         * Returns the row count.
         *
         * @return The row count.
         */
        public int getRowCount() {
            return this.data.length;
        }

        /**
         * Returns the value at the specified cell in the table.
         *
         * @param row  the row index.
         * @param column  the column index.
         *
         * @return The value.
         */
        public Object getValueAt(int row, int column) {
            return this.data[row][column];
        }

        /**
         * Sets the value at the specified cell.
         *
         * @param value  the value.
         * @param row  the row index.
         * @param column  the column index.
         */
        public void setValueAt(Object value, int row, int column) {
            this.data[row][column] = value;
            fireTableDataChanged();
        }

        /**
         * Returns the column name.
         *
         * @param column  the column index.
         *
         * @return The column name.
         */
        public String getColumnName(int column) 
        {
        	if(column > col_names.length - 1)
        		return null;
        	
        	return col_names[column];
        }
    }
}