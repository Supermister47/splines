import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Vector;
import java.util.Random;
import org.jblas.*;
import static org.jblas.FloatMatrix.*;
import static org.jblas.MatrixFunctions.*;



int screenWidth = 1000;
int screenHeight = 700;

int pointXDragged;
int pointYDragged;
boolean draggingPoint = false;

int POINT_SIZE = 14;


LinkedList<CubicPolynomial> splines;

Random rand;
int totalPoints;
int[] pointsScreen;

int chosenLineColor;
color lineColors[] = {color(77,238,234), color(116,238,21), color(255,231,0),
                      color(240,0,255), color(254,0,0), color(175,61,255),
                      color(255,59,148), color(166,256,41), color(22,133,248)};
color pointColors[];




class CubicPolynomial {
  public
    float a;
    float b;
    float c;
    float d;
    
    int xi;
    
  CubicPolynomial(float ai, float bi, float ci, float di, int x) {
    a = ai;
    b = bi;
    c = ci;
    d = di;
    xi = x;
  }
  
  
  float evaluate(int x) {
    
    float res = a;
    
    int aux = x - xi;
    res += b*aux;
    
    aux *= x - xi;
    res += c*aux;
    
    aux *= x - xi;
    res += d*aux;

    return res;
  }
}
  
  

class Point {
  public float x;
  public float y;
    Point(float c1, float c2) {
      x = c1;
      y = c2;
    }
}



// ---------------- Matrix operations ---------------- //

FloatMatrix forwardSubstitutionTridiag(FloatMatrix A, FloatMatrix b) {
  int m = A.rows;
  FloatMatrix x = new FloatMatrix(m);
  
  for (int i=m-1; 0 <= i; i--) {
    float sum = 0;
    for (int j=i+1; j < min(i+3, m); j++) {
      sum += A.get(i,j) * x.get(j);
    }
    
    x.put(i,0, (b.get(i) - sum) / A.get(i,i));
  }
  
  return x;
}



// Gaussian elimination optimized for tridiagonal matrices
FloatMatrix elimGaussianaTridiag(FloatMatrix A, FloatMatrix b) {
  float epsilon = 1e-15; //<>//
  int m = A.rows;
  FloatMatrix x = new FloatMatrix(m);
  
  for (int i=0; i < m-1; i++) {
    

    // If the element below the diagonal is 0 (detected with the epsilon) is setted as such, otherwise it would remain
    // as a value very close to 0
    if (abs(A.get(i+1, i)) < epsilon) {
      A.put(i+1,i, 0);
      continue;
    }
    
    // The element of the diagonal cannot be 0 because A is strictly diagonally dominant
    float pivot = A.get(i+1,i) / A.get(i,i);  
    
    for (int k=i; k < min(i+3, m); k++) {
      float res = A.get(i+1,k) - pivot * A.get(i,k);
      
      if (res != 0 && abs(res) < epsilon) {
        A.put(i+1,k, 0);
      }
      else {
        A.put(i+1,k, res);
      }
    }
    
    float resB =  b.get(i+1,0) - pivot * b.get(i,0);
      
    if (resB != 0 && abs(resB) < epsilon) {
      b.put(i+1,0, 0);
    }
    else {
      b.put(i+1,0, resB);
    }
     
  }
  
  // Once the matrix is triangular, the system is solved by doing forward substitution
  x = forwardSubstitutionTridiag(A, b);
  return x;
}
 //<>//



void piecewiseInterpolation() {
  splines.clear();
  
  // This auxiliary vector simplifies the access to Points
  Vector<Point> pointsList = new Vector<Point>();
  
  for (int x=0; x < screenWidth; x++) {
    int y = pointsScreen[x];
    
    if (y != 0) {
      pointsList.add(new Point(x, y));  
    }
  }
  
  int n = pointsList.size()-1;
  FloatMatrix coefMatrix = FloatMatrix.eye(n);
  FloatMatrix vectorB = new FloatMatrix(n);
  
  
  
  float oneThird = 1.0/3;
  float twoThirds = 2.0/3;
  
  // coefMatrix will be filled with the neccessary coeficents to calculate each of the splines
  for (int i=1; i < n; i++) {
    float bi;
    
    Point pi_minus_1 = pointsList.get(i-1);
    Point pi = pointsList.get(i);
    Point pi_plus_1 = pointsList.get(i+1);
    
    
    float hi = pi.x - pi_minus_1.x;
    float hi_plus_1 = pi_plus_1.x - pi.x;
    coefMatrix.put(i, i-1, oneThird*hi);
    coefMatrix.put(i, i, twoThirds*(hi + hi_plus_1));
    

    // If it's not the last row, the element is added to the right of the diagonal 
    if (i != n-1) {
      coefMatrix.put(i, i+1, oneThird*hi_plus_1);

      bi = (pi_plus_1.y - pi.y) / hi_plus_1 - (pi.y - pi_minus_1.y) / hi;
    }
    
    else {
      bi = (pi_plus_1.y - pi.y) / hi_plus_1 - (pi.y - pi_minus_1.y) / hi;
    }
    

    vectorB.put(i, 0, bi);
  }
  
  
  
  // Once the matrix is filled, we solve the system to obtain the quadratic coefficient of each spline
  FloatMatrix solC = elimGaussianaTridiag(coefMatrix, vectorB);
  //FloatMatrix solC = Solve.solve(coefMatrix, vectorB);
  
  float[] solA = new float[n];
  float[] solB = new float[n];
  float[] solD = new float[n];
  
  
  // Now the rest of the coefficients are calculated
  for (int i=0; i < n; i++) {

    Point pi = pointsList.get(i);
    Point pi_plus_1 = pointsList.get(i+1);
    
    
    float hi_plus_1 = pi_plus_1.x - pi.x;
    
    solA[i] = pi.y;
    
    if (i != n-1) {
      solB[i] = (((float)pi_plus_1.y - pi.y) / hi_plus_1) - hi_plus_1 * (twoThirds*solC.get(i) + oneThird*solC.get(i+1));
      solD[i] = oneThird * (solC.get(i+1) - solC.get(i)) / hi_plus_1;
    }
    else {
      solB[i] = (float)(pi_plus_1.y - pi.y) / hi_plus_1 - (twoThirds*solC.get(i)*hi_plus_1);
      solD[i] = -oneThird * (solC.get(i) / hi_plus_1);
    }
    
    splines.addLast(new CubicPolynomial(solA[i], solB[i], solC.get(i), solD[i], (int)(pi.x)));
  }
  
  /*
  print("Coeficientes ai\n");
  printArray(solA);
  print("\nCoeficientes bi\n");
  printArray(solB);
  print("\nCoeficientes di\n");
  printArray(solD);
  */
  
}







// ---------------- setup and draw ---------------- //

void setPointColors() {
  for (int i=0; i < lineColors.length; i++) {
    
    
    int red = round(red(lineColors[i]) * 0.5);
    int green = round(green(lineColors[i]) * 0.5);
    int blue = round(blue(lineColors[i]) * 0.5);
    
    pointColors[i] = color(red, green, blue);
  }
}



void createBackground() {
  background(10);
  strokeWeight(1);
  stroke(50);
  
  // Draw thinnest vertical lines
  for (int x=0; x < screenWidth; x+=50) {
    line(x, 0, x, screenHeight);
  }
  
  // Draw thinnest horizontal lines
  for (int y=0; y < screenHeight; y+=50) {
    line(0, y, screenWidth, y);
  }
  
  
  stroke(80);
  strokeWeight(2);
  
  // Draw thickest vertical lines
  for (int x=0; x < screenWidth; x+=100) {
    line(x, 0, x, screenHeight);
  }
  
  // Draw thickest horizontal lines
  for (int y=0; y < screenHeight; y+=100) {
    line(0, y, screenWidth, y);
  }
}



void settings() {
  size(screenWidth, screenHeight);
  rand = new Random();
}



void setup() {
  PImage icon = loadImage("icon.png");
  surface.setIcon(icon);
  
  splines = new LinkedList<CubicPolynomial>();

  pointsScreen = new int[screenWidth];
  totalPoints = 0;
  
  chosenLineColor = rand.nextInt(lineColors.length);
  pointColors = new int[lineColors.length];
  setPointColors();
  stroke(lineColors[chosenLineColor]); //<>//
}



void draw() {
  
  createBackground();
  
  // Set line and point colors
  stroke(lineColors[chosenLineColor]);
  fill(pointColors[chosenLineColor], 255f);
  
  int countPoints = 0;
  int lastXPoint = -1;
  
  // The first point (left to right) is obtained
  int firstPoint_x = -1;
  for (int x=0; x < screenWidth; x++) {
    int y = pointsScreen[x];
    
    if (y != 0) {
      firstPoint_x = x;
      lastXPoint = x;
      countPoints++;

      break;
    }
  }
  
  
  // If there is more than one point, the interpolation is computed and drawed on the screen
  if (totalPoints > 1) {
    piecewiseInterpolation();
    
    
    int x = firstPoint_x;
    while (x < screenWidth && countPoints < totalPoints) {
      
      lastXPoint = x;
      
      int y = pointsScreen[x];
      do {    
        int thisPointValue = round(splines.get(countPoints-1).evaluate(x));
        int nextPointValue = round(splines.get(countPoints-1).evaluate(x+1));
        line(x, thisPointValue, x+1, nextPointValue);
        
        x++;
        y = pointsScreen[x];
      } while (y == 0);
      
      ellipse(lastXPoint, pointsScreen[lastXPoint], POINT_SIZE, POINT_SIZE);
      countPoints++;
    }
    
    ellipse(x, pointsScreen[x], POINT_SIZE, POINT_SIZE);
  }
  
  // Otherwise, the only point is shown
  else if (totalPoints == 1) {
      ellipse(firstPoint_x, pointsScreen[firstPoint_x], POINT_SIZE, POINT_SIZE);
  }
  
  // If the user has selected a Point, the coordinates of the position are displayed
  if (draggingPoint) {
    String coordinateText = "(";
    coordinateText += (float)pointXDragged/10 + ", ";
    coordinateText += nf((float)screenHeight/10-(float)pointsScreen[pointXDragged]/10, 0, 1) + ")";
    textSize(20);
    fill(250);
    text(coordinateText, 50, 50);
  }
}




// ---------------- Input ---------------- //

void keyReleased() {
    if  (key == DELETE && draggingPoint) {
      pointsScreen[pointXDragged] = 0;
      pointXDragged = -1;
      totalPoints--;
      draggingPoint = false;
    }
    
    // Change colors
    else if (key == CODED) {
      if (keyCode == RIGHT) {
        chosenLineColor = Math.floorMod(chosenLineColor+1, lineColors.length);
      }
      else if (keyCode == LEFT) {
        chosenLineColor = Math.floorMod(chosenLineColor-1, lineColors.length);
      }
    }
}

void mouseDragged() {
  
  int mouseXConstrained = constrain(mouseX, 0, screenWidth-1);
  int mouseYConstrained = constrain(mouseY, 1, screenHeight-1);

  if (draggingPoint && (pointsScreen[mouseXConstrained] == 0 || mouseXConstrained == pointXDragged)) {
    pointsScreen[pointXDragged] = 0;
    pointsScreen[mouseXConstrained] = mouseYConstrained;
    pointXDragged = mouseXConstrained;
  }
}



void mousePressed() {
  
  draggingPoint = true;
  //<>//
  int mouseXConstrained = constrain(mouseX, 0, screenWidth-1);
  int mouseYConstrained = constrain(mouseY, 1, screenHeight-1);
  

  // Select a nearby point, if there is one close to the mouse
  int pointXd;
  for (int d=0; d < POINT_SIZE; d++) {
    
    pointXd = mouseXConstrained-d;
    if (pointsScreen[max(pointXd, 0)] != 0) {
      
      if (pointsScreen[pointXd] - POINT_SIZE <= mouseYConstrained && mouseYConstrained <= pointsScreen[pointXd] + POINT_SIZE) {
        mouseYConstrained = pointsScreen[pointXd];
        mouseXConstrained = pointXd;
        pointsScreen[pointXd] = 0;
        totalPoints--;
        break;
      }  
    }

    else if (pointsScreen[min(mouseXConstrained+d, screenWidth-1)] != 0) {
      pointXd = mouseXConstrained+d;
      if (pointsScreen[pointXd] - POINT_SIZE <= mouseYConstrained && mouseYConstrained <= pointsScreen[pointXd] + POINT_SIZE) {
        mouseYConstrained = pointsScreen[pointXd];
        mouseXConstrained = pointXd;
        pointsScreen[pointXd] = 0;
        totalPoints--;
        break;
      }  
    }
  }
  
  // If the mouse it is not near a point but it is in the same x coordinate as some point, then
  // no point is added
  if (pointsScreen[mouseXConstrained] != 0) {
    return;
  }
  
  pointsScreen[mouseXConstrained] = mouseYConstrained;
  pointXDragged = mouseXConstrained;
  pointYDragged = mouseYConstrained;
  
  totalPoints++;

  return;
}



void mouseReleased() {
  draggingPoint = false;
}
