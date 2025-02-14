import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLineEdit, QLabel
from PyQt5.QtCore import Qt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np
from units import *




# Your existing function to generate the interaction diagram
def generate_interaction_diagram(b,h,CC_to_CL_long,fpc, num_lay, bar_per_lay, bar_size, fy):
    
    Ec = 57000*(fpc/psi)**0.5*psi
    εcu = 0.003
    ß1 = 0.85-(0.05*(fpc-4000*psi))/1000
    
    def create_rebar_array(num_lay,bar_per_lay,bar_size,CC_to_CL_long):
        rebar_array = []
        for i in range(num_lay):
            rebar_array.append([bar_per_lay,bar_size,CC_to_CL_long+i*(h-2*CC_to_CL_long)/(num_lay-1)])
        
        return rebar_array

    def find_strain(Z,rebar_array,εcu):
        Es = 29000*ksi
        εs1 = Z*fy/Es
        c = c = εcu*((h-CC_to_CL_long)/(εcu-εs1))
        a = c*ß1
        if a > h:
            a = h
            
        #print(εs1)
        
        if εs1 <= -εcu-fy/Es:
            φ = 0.9
        elif εs1 < -fy/Es and εs1 > -fy/Es-εcu:
            φ = 0.65+0.25*(-εs1-fy/Es)/εcu
        elif εs1 < 0 and εs1 >= -fy/Es:
            φ = 0.65
        else:
            φ = 0.65


        strains = []

        for reb in rebar_array:

            strains.append(εcu*(c-reb[2])/c)

        return strains,a,φ

    def find_rebar_force(strains,rebar_array,fpc=4*ksi,fy=60*ksi,Es=29000*ksi):
        forces = []
        #print(len(strains))
        for strain_count in range(len(strains)):

            if strains[strain_count] >= fy/Es:
                forces.append((fy-0.85*fpc)*rebar_array[strain_count][0]*rebars[rebar_array[strain_count][1]]['area'])

            elif strains[strain_count] <= -fy/Es:
                forces.append((-fy)*rebar_array[strain_count][0]*rebars[rebar_array[strain_count][1]]['area'])

            elif strains[strain_count] >0:
                forces.append(((strains[strain_count]*Es)-0.85*fpc)*rebar_array[strain_count][0]*rebars[rebar_array[strain_count][1]]['area'])

            else:
                forces.append((strains[strain_count]*Es*rebar_array[strain_count][0])*rebars[rebar_array[strain_count][1]]['area'])

        return forces

    def find_PM(a,b,h,forces, rebar_array,φ, fpc=4*ksi):
        """
        This generates a singular point for PM diagram in KIP, KIP-FT
        """
    
        C_c = 0.85*fpc*a*b #Compressive force Kips

        Pn = sum(forces)+C_c
        
        Mn = []
        for force_count in range(len(forces)):
            Mn.append(forces[force_count]*(h/2-rebar_array[force_count][2]))
                        
        Mn = sum(Mn)+C_c*((h/2)-a/2)
        #return Pn*φ,Mn*φ/12
        return Pn*φ,Mn*φ/12

    def find_all_Z_PM(rebar_array,εcu,b,h,fpc,fy=4*ksi):
        PM_points = []
        for Z in np.arange(1,-1000.1,-0.1):
            strains,a,φ = find_strain(Z,rebar_array,εcu)
            forces = find_rebar_force(strains,rebar_array)
            PM_points.append(find_PM(a,b,h,forces,rebar_array,φ))
            
        return PM_points
    
    def plot_truncated_PM(PM_points,rebar_array, b, h, fpc, fy):
        As = 0
        for reb in rebar_array:
            As += reb[0]*rebars[reb[1]]['area']
            
        Max_φPn = 0.8*0.65*(0.85*fpc*((b*h)-(As))+(As*fy))
        print(Max_φPn)
        #PM_points = np.transpose(PM_points)
        #PM_points[0] = np.where(PM_points[0]<Max_φPn,PM_points[0],Max_φPn)
        PM_trunc = np.where(PM_points[0]<Max_φPn,PM_points[0],Max_φPn)
        x = PM_points[1]
        y = PM_trunc

        #fig, ax = plt.subplots()
        ##ax.plot(PM_points[1],PM_points[0])
        #ax.plot(PM_points[1],PM_trunc)
        #ax.grid()
        #ax.set(xlabel='Moment (Kip-ft)', ylabel='Axial Loading (Kip)',
        #title='ACI 318 - PM Interaction Diagram')
        #plt.show()  
        return x, y

    rebar_array = create_rebar_array(num_lay, bar_per_lay, bar_size,CC_to_CL_long)
    x,y = plot_truncated_PM(np.transpose(find_all_Z_PM(rebar_array,εcu,b,h,fpc,fy=4*ksi)),rebar_array, b, h, fpc, fy)
    
    return x, y


class MyWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Reinforced Concrete Column Interaction Diagram")
        self.setGeometry(100, 100, 800, 600)

        # Set up the layout
        self.layout = QVBoxLayout()

        # Input fields for geometry (width, depth)
        self.geometry_layout = QHBoxLayout()
        self.geometry_label = QLabel("Enter Geometry: b (in), h (in), CC to Long (in)")
        self.width_input = QLineEdit()
        self.depth_input = QLineEdit()
        self.CC_to_CL_long_input = QLineEdit()
        self.geometry_layout.addWidget(self.geometry_label)
        self.geometry_layout.addWidget(self.width_input)
        self.geometry_layout.addWidget(self.depth_input)
        self.geometry_layout.addWidget(self.CC_to_CL_long_input)
        self.layout.addLayout(self.geometry_layout)

        # Input fields for reinforcement (fpc, fy)
        self.reinforcement_layout = QHBoxLayout()
        self.reinforcement_label = QLabel("Enter Reinforcement: num_lay, bar_per_lay, bar_size, fpc (psi), fy (ksi)")
        self.num_lay_input = QLineEdit()
        self.bar_per_lay_input = QLineEdit()
        self.bar_size_input = QLineEdit()
        self.fpc_input = QLineEdit()
        self.fy_input = QLineEdit()
        self.reinforcement_layout.addWidget(self.reinforcement_label)
        self.reinforcement_layout.addWidget(self.num_lay_input)
        self.reinforcement_layout.addWidget(self.bar_per_lay_input)
        self.reinforcement_layout.addWidget(self.bar_size_input)
        self.reinforcement_layout.addWidget(self.fpc_input)
        self.reinforcement_layout.addWidget(self.fy_input)
        self.layout.addLayout(self.reinforcement_layout)

        # Create a button that will generate the plot
        self.button = QPushButton("Generate Diagram", self)
        self.button.clicked.connect(self.generate_and_plot)
        self.layout.addWidget(self.button)

        # Create a Matplotlib Figure and Canvas
        self.figure = plt.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)

        self.ax = self.figure.add_subplot(111)  # Create subplot once



        self.setLayout(self.layout)

    def generate_and_plot(self):
        # Get user input values
        try:
            b = float(self.width_input.text())*inch
            h = float(self.depth_input.text())*inch
            CC_to_CL_long = float(self.CC_to_CL_long_input.text())*inch
            fpc = float(self.fpc_input.text())*psi
            num_lay = int(self.num_lay_input.text())
            bar_per_lay = int(self.bar_per_lay_input.text())
            bar_size = str(self.bar_size_input.text())
            fy = float(self.fy_input.text())*ksi

            
        except ValueError:
            print("Invalid input values!")
            return

        # Generate the interaction diagram using the existing function
        #geometry = (b, h)
        #reinforcement = (fpc, fy)
        #x, y = generate_interaction_diagram(geometry, reinforcement)
        x, y = generate_interaction_diagram(b,h,CC_to_CL_long,fpc, num_lay, bar_per_lay, bar_size, fy)

        # Plot the generated diagram

    
        #ax = self.figure.add_subplot(111)
        self.ax.clear()
        self.ax.plot(x, y)
        self.ax.relim()
        self.ax.autoscale_view()
        self.ax.grid()
        self.ax.set(xlabel='Moment (Kip-ft)', ylabel='Axial Loading (Kip)',title='ACI 318 - PM Interaction Diagram')
        self.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())

