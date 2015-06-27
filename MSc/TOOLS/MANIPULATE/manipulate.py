#!/usr/bin/env python

import sys
import os
import subprocess
import pygtk
import gtk
import gtk.glade
import time


#### Clases ###
class gnuplot:
    """Simple pipe to a gnuplot process"""
    command="gnuplot"
    def __init__(self):
        self.gnuplot = os.popen(self.command, 'w')
        #self.process = subprocess.Popen(self.command,
        #                                shell=True,
        #                                stdin=subprocess.PIPE)
        #self.gnuplot = self.process.stdin
        #self.out = self.process.stdout
        #self.err = self.process.stderr
        self.write = self.gnuplot.write
        self.flush = self.gnuplot.flush

    def close(self):
        if self.gnuplot is not None:
            self.gnuplot.close()
            self.gnuplot = None

    def __del__(self):
        self.close()

    def __call__(self, s):
        """Send a command string to gnuplot, followed by newline."""
        self.write(s + '\n')
        self.flush()


class ManipulateGTK:
    """The main window of the manipulate process"""
    def __init__(self, filenames):
        #Set the Glade file
        home=os.environ.get('HOME')
        gladefile=home+"/RESEARCH/PROGRAMMING/TOOLS/MANIPULATE/manipulate.glade"  
        self.wTree = gtk.glade.XML(gladefile)
        self.index=0

        self.plots=[]
        for file in filenames:
            self.plots.append( gnuplot() )

        # Start the gnuplot connections
        for i in range(len(filenames)):
            self.plots[i]("set terminal wxt")
            self.plots[i]("_i="+str(self.index))
            file=open(filenames[i],"r")
            for line in file:
                self.plots[i](line)
            file.close()

        #Create our dictionay and connect it
        dic = { "on_window1_destroy" : gtk.main_quit ,
                "on_Plot_clicked" : self.Plot, 
                "on_Play_clicked" : self.Play,
                "on_Refresh_clicked" : self.Refresh,
                "on_Barra_value_changed" : self.Barra
              }
        self.wTree.signal_autoconnect(dic)
    
    def Plot(self, widget):
        for gpipe in self.plots:
            gpipe("_i="+str(self.index))
            gpipe("replot")
            #gpipe("show xrange")
            #     print gpipe.err

    def Play(self,widget):
        for i in range(self.wTree.get_widget('Barra').get_property("adjustment").upper + 1):
            self.index=i
            self.Plot(widget)
            #time.sleep(0.5)

    def Barra(self,widget):
        self.index=round(self.wTree.get_widget('Barra').get_property("adjustment").value,0)
        self.Plot(widget)

    def Refresh(self, widget):
        self.wTree.get_widget('Barra').get_property("adjustment").upper=\
            self.wTree.get_widget('upper').get_property("adjustment").value
        self.wTree.get_widget('Barra').get_property("adjustment").lower=\
            self.wTree.get_widget('lower').get_property("adjustment").value



if __name__ == "__main__":
    a=ManipulateGTK(sys.argv[1:])
    gtk.main()

