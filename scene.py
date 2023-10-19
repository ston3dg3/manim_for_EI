from manim import *
import scipy.integrate as integrate
import numpy as np



class CreateCircle(Scene):
    def construct(self):
        circle = Circle() #creates a circle
        circle.set_fill(PINK, opacity=0.5) # set color and transparency
        self.play(Create(circle)) #show circle on screen





def normalize(v): return v / (np.linalg.norm(v))
def norm(u, v): return np.linalg.norm(np.dot(u,v))

# Projektion von Vektor v entalang (basis)-Vektor u (u ungleich 0)
def proj(v, u):
    if(np.all(u==0)):
        # TODO: raise exception: Basisvektor ist Null
        return
    else:
        return (np.dot(v,u))*u
        # Die Vektoren u und [v-proj(v,u)] sind orthogonal


# Ziel: Gram Schmidt für alle s_i 

# V ist eine menge von m Vektoren mit beliebiger Basis (v_1,v_2,...,v_m) 
# Nach Gram Schmidt, bilden die Vektoren U:=(u_1,u_2,...,u_k) mit k >= m, eine Orthonormale Basis
def GramSchmidt(V):
    
    # Alle Nullvektoren entfernen - damit stellen wir sicher dass V[0] immer ein echtes Vektor ist
    V = V[~np.all(V == 0, axis=1)]

    # Anzahl von Nicht-Null Vektoren in V
    M = len(V)
    # Basis U von V erstellen mit N Vektoren (maximal N <= M Vektoren)
    # wenn M = N, dann ist V schon eine Orthogonale Basis!!!
    U = []
    
    for v in V:
        v_u = []
        for u in U:
            v_u.append(proj(v,u))
            # Draw animation Vector
        # hier haben wir die Projektionen v_u von v auf alle bisherige basisvektoren u von U 
        u = v - sum(v_u)
        if (u > 1e-10).any():  
            U.append(normalize(u))

    
        
    U_basis = np.array(U)
    return U_basis


V =np.array([   [3,0,4],
                [1,0,0],
                [1,0,2]])

U_basis = GramSchmidt(V)
print (U_basis)




class Vectors(VectorScene):
    def construct(self):
        
        plane = self.add_plane(animate = True).add_coordinates()
        vector = self.add_vector([-3,-2], color = YELLOW)

        basis = self.get_basis_vectors()
        self.add(basis)
        self.vector_to_coords(vector=vector)

        vector2 = self.add_vector([2,2])
        self.write_vector_coordinates(vector=vector2)

        



class Axes3DExample(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()
        self.add_plane(animate = True)


        x_label = axes.get_x_axis_label(Tex("x"))
        y_label = axes.get_y_axis_label(Tex("y")).shift(UP * 1.8)

        # 3D variant of the Dot() object
        dot = Dot3D()

        # zoom out so we see the axes
        self.set_camera_orientation(zoom=0.5)

        self.play(FadeIn(axes), FadeIn(dot), FadeIn(x_label), FadeIn(y_label))

        self.wait(0.5)

        # animate the move of the camera to properly see the axes
        self.move_camera(phi=75 * DEGREES, theta=30 * DEGREES, zoom=1, run_time=1.5)

        # built-in updater which begins camera rotation
        self.begin_ambient_camera_rotation(rate=0.15)

        # one dot for each direction
        upDot = dot.copy().set_color(RED)
        rightDot = dot.copy().set_color(BLUE)
        outDot = dot.copy().set_color(GREEN)

        self.wait(1)

        self.play(
            upDot.animate.shift(UP),
            rightDot.animate.shift(RIGHT),
            outDot.animate.shift(OUT),
        )

        self.wait(2)

class ThreeDSceneExample(ThreeDScene):
    def construct(self):
        axes = ThreeDAxes()  # Create 3D axes
        self.set_camera_orientation(phi=60 * DEGREES, theta=-45 * DEGREES)

        self.play(Create(axes))
        self.begin_ambient_camera_rotation(0.1)  # Spin camera

        # Define vectors

        # V =np.array([   [3,0,4],
        #                 [1,0,0],
        #                 [1,0,2]])
        
        V =np.array([   [3,0,5],
                        [1,0,-2],
                        [-1.5,3,2]])

        V = V[~np.all(V == 0, axis=1)]
        U = []
        x = 0

        for v in V:
            x = x + 1
            vector_text = MathTex("\\vec{v}_"+str(x), color=WHITE).move_to(v)
            self.play(Create(Arrow3D(ORIGIN, v, color=YELLOW)), Write(vector_text))
            v_u = []
            y=0

            for u in U:
                pr = proj(v,u)
                v_u.append(pr)
                # Draw animation Vector
                y = y + 1
                strr = f"{x}{y}"
                projVector_text = MathTex("\\vec{v}_{"+strr+"}", color=WHITE).move_to(pr * 0.5)
                self.play(Create(Arrow3D(ORIGIN, pr, color=random_bright_color())), Write(projVector_text))

            # hier haben wir die Projektion v_u von v auf alle bisherige basisvektoren u von U 
            u = v - sum(v_u)
            if (u > 1e-10).any():  
                U.append(normalize(u))
                # Draw New Basis vector
                vector_text = MathTex("\\vec{u}_"+str(x), color=WHITE).move_to(u)
                self.play(Create(Arrow3D(ORIGIN, u, color=RED)), Write(vector_text))
                self.wait(1.5)
            

        self.wait(5)
        self.stop_ambient_camera_rotation()
        self.wait(2)