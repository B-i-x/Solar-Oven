
import math as m
from re import X
import numpy as np
import plotly.express as px

from plotly.express import data



class Ovens:
    '''this class creates objects of ovens
    ovens have the following variables:

        completely static variables (variables you cannot change):
            aUw is the x^2 coeff of the trendline for the mylar sheet
            bUw is the x coeff of the trendline for the mylar sheet 
            cUw is the constant coeff of the trendline for the mylar sheet

        input variables:
            tao is the optical transmission coefficient of the window material
            a is the absorption factor of the cavity walls (0 ≤ a ≤ 1)
            Tambient is the ambient temperature in celsius
            SolarPower is the solar power
            n is the number of mylar sheets, only can be 1 or 2
            moverw is the ratio of the length of the panels to the width of the box, 0 < x < 4
            W is width of the box in meters 
            S is the space between the inner box and the other box in meters
            G is the energy gain of the reflectors
            r is the reflectiveness of the aluminum foil


    '''

    def __init__(self):
        self.Tambient = 26 #degress celsius
        self.SolarPower = 1000
        self.angle = 90
        self.n = 2
        self.r = 0.70
        self.moverw = 3
        self.tao = .92
        self.a = 0.90

        self.W = "default"
        self.H = "default"
        self.S = "default"

        self.G = None

    def setdimensions(self):

        if self.W == "default":
            self.W = .1
        
        if self.H == "default":
            self.H = .001/(self.W**2)

        if self.S == "default":
            self.S = .05


    def innerbox(self):
        self.WindowArea = (self.W)**2
        self.Inner5Area = 4*self.W*self.H + self.WindowArea

    def HeatTransferm(self):
        InnerCardBoardWall = 1
        ICBWThickness = InnerCardBoardWall *.004
        OuterCardBoardWall = 1
        OCBWThickness = OuterCardBoardWall * .004
        InsulationThickness = self.S
        # k values are the thermal conductivity of the materials
        Cardboardk = .064
        Newspaperk = .123
        self.HeatTransfer = 1/((ICBWThickness/Cardboardk) + (InsulationThickness/Newspaperk) + (OCBWThickness/Cardboardk))

    def gmath(self):
        moverl = self.moverw
        #print(moverl)
        alpha = m.degrees(m.asin((-moverl+(((moverl**2)+8)**0.5))/4))
        self.alpha = alpha
        #alpha = 16.31 # degrees
        NumReflectors = 4
        #pabsorbed = SolarPower * WindowArea * a * (tao**n)
        if self.G == None:
            #self.G = (1+NumReflectors*self.r*moverl*m.sin(m.radians(alpha)))*(1+moverl*m.sin(m.radians(alpha)))
            #self.G = 1+NumReflectors*self.r*(moverl)*m.sin(m.radians(alpha))+NumReflectors*self.r*((moverl*self.W)**2)*(m.sin(m.radians(alpha)))**2
            self.G = 1 + NumReflectors*self.r*moverl*m.sin(m.radians(alpha))
            
    def temp(self):
        if self.n == 1:
            aUw = 0.0008005894
            bUw = .00050596
            cUw = 6.5899236
        if self.n == 2 :
            aUw = 0.0003925812
            bUw = float(-0.0015663384)
            cUw = 3.2775005875
        #converting to form ax^3+bx^2+cx+d=0
        # a lvoneq = a in the form of the cubic formula and since it is simpler with two levels it is level one, https://math.vanderbilt.edu/schectex/courses/cubic/)
        alvone = aUw*self.WindowArea 
        blvone = bUw*self.WindowArea - self.Tambient*aUw*self.WindowArea
        clvone = (cUw*self.WindowArea + self.HeatTransfer*self.Inner5Area)-(self.Tambient*bUw*self.WindowArea)
        dlvone = -(self.G*self.SolarPower*self.WindowArea*(self.tao**self.n)*self.a*m.sin(m.radians(self.angle)))-(self.Tambient*self.HeatTransfer*self.Inner5Area)-(self.Tambient*cUw*self.WindowArea)

        #using numpy
        coeff = [alvone, blvone, clvone, dlvone]
        #print(self.coeff)
        roots = np.roots(coeff)
        #print(roots)
        #print(alvone, blvone, clvone, dlvone)

        # printing the real root
        realindex = np.where(np.isreal(roots))
        self.Tinternal = round(float((roots[realindex]).real), 5)
        return self.Tinternal

    def cost(self):
        extra = .15
        outsidebox = 4*((self.W+self.S*2)*(self.S+self.H))+(self.W+self.S*2)**2+((self.W+self.S*2)**2-self.WindowArea)
        #print(outsidebox)
        #reflector area math
        phi = 90-self.alpha
        M = self.moverw*self.W

        diagonal = (((2**(0.5)*M*m.cos(m.radians(phi)))**2)+(M*m.sin(m.radians(phi)))**2)**(0.5)
        #print(diagonal)
        #print(M)
        extraonsides = (diagonal**2-M**2)**(0.5)
        #print(extraonsides)
        top = (2*extraonsides)+self.W
        #print(top)
        reflectorarea = 4*((extraonsides+self.W)*M)
        #print(reflectorarea)

        cardboardcost = 1.75*((outsidebox+self.Inner5Area+reflectorarea)+((outsidebox+self.Inner5Area+reflectorarea)*extra))
        #print(cardboardcost)
        mylarcost = self.n*.25
        aluminumfoilcost = reflectorarea*.55
        totalcost = round(cardboardcost+mylarcost+aluminumfoilcost, 2)

        self.tcost = totalcost
        return totalcost

    
    def create(self):
        self.setdimensions()
        self.innerbox()
        self.gmath()
        self.HeatTransferm()
        self.temp()
        self.cost()

        self.tempperdollar = round(self.Tinternal/self.tcost, 5)



class OvenTest:
    '''this class tests the ovens in whatever way you want to'''
    def __init__(self, steps, x=None, y=None, z=None) -> None:
        '''this is the input the test case you want on the graph, must be at least two varaiables with their ranges'''

        self.steps = steps

        self.inputs = 0 #of variables

        if x is not None:
            self.x_var_name = x[0]
            self.x_input_list = self.create_input_list(x[1], x[2])
            self.inputs += 1

        if y is not None:
            self.y_var_name = y[0]
            self.y_input_list = self.create_input_list(y[1], y[2])
            self.inputs += 1

        if z is not None:
            self.z_var_name = z[0]
            self.z_input_list = self.create_input_list(z[1], z[2])
            self.inputs += 1

    def create_input_list(self, min, max):

        nump_array = np.linspace(min, max, self.steps)

        return list(nump_array)

    def labels(self, labels):

        self.label_list = labels

    def file_name(self, filename):

        self.filename = filename

    def set_constants(self, setting=None, varlist=None):
        
        if setting == "Default":
            pass

        self.varlist = varlist

    def ouput(self, result, result2=None):

        self.r1 = result
        self.r2 = result2

    def one_oven(self):

        single_instance = Ovens()

        for var_values in self.varlist:
            setattr(single_instance, var_values[0], var_values[1])

        single_instance.create()

        result_single = getattr(single_instance, self.r1)

        print(vars(single_instance))

        return result_single

    def oven_creation(self):

        x_output = []
        y_output = []
        z_output = []

        result_list = []
        result2_list = []


        if self.inputs == 1:
            permutation = [[x] for x in self.x_input_list]
        elif self.inputs == 2:
            permutation = [[x,y] for x in self.x_input_list for y in self.y_input_list]
        elif self.inputs == 3:
            permutation = [[x,y,z] for x in self.x_input_list for y in self.y_input_list for z in self.z_input_list]

        self.instances = [Ovens() for i in range(self.steps**self.inputs)]

        for index, instance in enumerate(self.instances):

            if self.varlist is not None:
                for var_values in self.varlist:
                    setattr(instance, var_values[0], var_values[1])

            x = permutation[index][0]
            setattr(instance, self.x_var_name, x)
            x_output.append(x)

            if self.inputs > 1:
                
                y = permutation[index][1]
                setattr(instance, self.y_var_name, y)
                y_output.append(y)

                if self.inputs > 2:
                    
                    z = permutation[index][2]
                    setattr(instance, self.z_var_name, z)
                    z_output.append(z)

            instance.create()

            result_list.append(getattr(instance, self.r1))

            if self.r2 is not None:
                result2_list.append(getattr(instance,self.r2))

        data = [x_output, y_output, z_output, result_list, result2_list]

        for dataset in data:
            if len(dataset) == 0:
                data.remove(dataset)

        return data


    def graph(self, df, data, labels):

        for dataset in data:
            if len(dataset) == 0:
                data.remove(dataset)

        datasets = len(data)


        if datasets == 2:
            fig = px.scatter(df, x=data[0], y=data[1], 
            labels={
                "x" : labels[0],
                "y" : labels[1]
            },
            title=f"{labels[0]} and its effect on {labels[1]}" )

        elif datasets == 3:
            fig = px.scatter_3d(df, x=data[0], y=data[1], z=data[2],
            labels={
                "x" : labels[0],
                "y" : labels[1],
                "z" : labels[2]
            },
            title=f"{labels[0]} and {labels[1]} its effect on {labels[2]}" )
        
        elif datasets == 4:
            fig = px.scatter_3d(df, x=data[0], y=data[1], z=data[2], color=data[3],
            labels={
                "x" : labels[0],
                "y" : labels[1],
                "z" : labels[2],
                "color" : labels[3]
            },
            title=f"{labels[0]} and {labels[1]} and {labels[2]} and its effect on {labels[3]}" )

        fig.show()

        filepath = r"C:\Users\alexr\Documents\Solar-Oven-Eng102-Python-Code\Solar-Oven-Eng102-Python-Code\{}.html".format(self.filename)

        fig.write_html(filepath)


    def print_solution(self, data=None, labels=None, singlevalue=None):
        
        if self.inputs != 0:

            for dataset in data:

                dataset_index = data.index(dataset)

                max_index = dataset.index(max(dataset))

                print(f"the max of {labels[dataset_index]} is {dataset[max_index]}")

                print(f"the settings are \n{vars(self.instances[max_index])}")

        else:
            print(f"the {self.r1} is {singlevalue}")

    def run_test(self):

        if self.inputs != 0:

        
            data = self.oven_creation()
            
            df = px.data.iris()

            self.graph(df=df, data=data, labels=self.label_list)

            self.print_solution(data=data, labels=self.label_list)

        else:
            singe_value = self.one_oven()

            self.print_solution(singlevalue=singe_value)




def widthandspace():
    '''this goes through width and space at the same time and graphs the resultant temperature and temp/dollar'''
    x = ["W", .01, .15]
    y = ["S", .01, .15]
    data_specifity = 25

    set_variables = [
        ["G", 4]
    ] 
    testcase = OvenTest(data_specifity, x, y)
    
    testcase.set_constants(setting="Default", varlist=set_variables)

    testcase.labels(labels=["Width of Box", "Space Between Inner and Outer", "Temperature"])
    testcase.ouput("Tinternal")
    testcase.file_name("width, space, temperature")
    testcase.run_test()

def static1():
    '''static'''
    variables = [
        ["W", .1005050505050504],
        ["H", .001/(.1005050505050504**2)],
        ["S", .04818181818181818],
        ["Tambient", 26],
        ["SolarPower", 1000],
        ["moverw", 3],
        ["G", 2]
    ]
    test2case = OvenTest(0)
    test2case.set_constants(varlist=variables)
    test2case.ouput("Temp/cost")
    test2case.run_test()

def wandtempcost():
    x = ["W", .01, .15]
    data_specifity = 100
    testcase = OvenTest(data_specifity, x)
    testcase.ouput("Temp/cost")
    testcase.set_constants(setting="Default")
    testcase.run_test()

def wandtemp():
    x = ["W", .01, .15]
    data_specifity = 200
    testcase = OvenTest(data_specifity, x)

    testcase.labels(["Width of Box", "Temperature"])
    testcase.ouput("Tinternal")
    testcase.set_constants(setting="Default")
    testcase.file_name("Width, temp")

    testcase.run_test()

def threeinputs():
    x = ["W", .01, .15]
    y = ["S", .01, .15]
    z = ["G", 2, 8]
    data_specifity = 25

    testcase = OvenTest(data_specifity, x, y, z)
    
    testcase.set_constants(setting="Default")

    testcase.labels(labels=["Width of Box", "Space Between Inner and Outer", "G value", "Temperature"])
    testcase.ouput("Tinternal")
    testcase.file_name("width, space, G, temperature")
    testcase.run_test()


widthandspace()
#static1()
#wandtempcost()
#wandtemp()
#threeinputs()



