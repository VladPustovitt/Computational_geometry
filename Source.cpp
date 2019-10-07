#include<iostream>
#include<cmath>
#include<iomanip>
#include<algorithm>
#include<vector>

using namespace std;
const double PI = acos(-1);

class Point {
private:
	double x, y;
public:
	Point(double x = 0, double y = 0) {
		this->x = x;
		this->y = y;
	}
	Point(const Point& p) {
		x = p.x;
		y = p.y;
	}
	Point& operator= (const Point& p) {
		x = p.x;
		y = p.y;
		return *this;
	}
	double findPolarFi() {
		double result = atan2(y, x);
		if (result < 0)
			result += 2 * PI;
		return result;
	}
	double findDistance(Point p) {
		return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
	}
	bool operator == (const Point& p) {
		return (x == p.x) && (y == p.y);
	}
	double X() {
		return x;
	}
	double Y() {
		return y;
	}
	Point operator - () {
		return Point(-x, -y);
	}
	friend istream& operator >> (istream&, Point&);
	friend ostream& operator << (ostream&, const Point&);
	friend class Line;
	friend class Circle;
	friend class Vector;
};
istream& operator >> (istream& in, Point& p) {
	in >> p.x >> p.y;
	return in;
}
ostream& operator << (ostream& out, const Point& p) {
	out << p.x << ' ' << p.y;
	return out;
}

class Vector {
private:
	double x, y;
public:
	Vector(double x = 0, double y = 0) {
		this->x = x;
		this->y = y;
	}
	Vector(const Vector& v) {
		x = v.x;
		y = v.y;
	}
	Vector(const Point& p1, const Point& p2) {
		x = p2.x - p1.x;
		y = p2.y - p1.y;
	}
	Vector& operator=(const Vector& v) {
		x = v.x;
		y = v.y;
		return *this;
	}
	double size() {
		return sqrt(x * x + y * y);
	}
	Vector operator * (double a) {
		return Vector(x * a, y * a);
	}
	Vector operator / (double a) {
		return Vector(x / a, y / a);
	}
	Vector operator + (Vector v) {
		return Vector(x + v.x, y + v.y);
	}
	Vector operator - (Vector v) {
		return Vector(x - v.x, y - v.y);
	}
	Vector operator - () {
		return Vector(-x, -y);
	}
	double operator () (Vector v) {
		return x * v.x + y * v.y;
	}
	double operator [] (Vector v) {
		return x * v.y - y * v.x;
	}
	friend istream& operator >> (istream&, Vector&);
	friend ostream& operator << (ostream&, const Vector&);
	friend class Line;
	friend class Circle;
	friend class Point;
	friend Point& ParallelCarry(const Point&, const Vector&);
};
double findFi(Vector v1, Vector v2) {
	double result = acos((v1(v2) / (v1.size() * v2.size())));
	return result;
}
istream& operator>>(istream& in, Vector& v){
	in >> v.x >> v.y;
	return in;
}
ostream& operator<<(ostream& out, const Vector& v) {
	out << v.x << ' ' << v.y;
	return out;
}

class Line {
private:
	double a, b, c;
public:
	Line(double a = 1, double b = 0, double c = 0) {
		this->a = a;
		this->b = b;
		this->c = c;
	}
	Line(const Point& p1, const Point& p2) {
		Vector p1p2(p1, p2);
		a = p1p2.y;
		b = -p1p2.x;
		c = -a * p1.x - b * p1.y;
	}
	Line(const Point& p, const Vector& v) {
		a = v.x;
		b = v.y;
		c = -a * p.x - b * p.y;
	}
	Line(const Line& A) {
		a = A.a;
		b = A.b;
		c = A.c;
	}
	Line& operator=(const Line& A) {
		a = A.a;
		b = A.b;
		c = A.c;
		return *this;
	}
	bool isPointOnLine(const Point& p) {
		return  a * p.x + b * p.y + c == 0;
	}
	short int isnotPointOnLine(const Point& p) {
		return a * p.x + b * p.y + c > 0 ? 1 : -1;
	}
	double PointDistance(const Point& p) {
		double D = abs(a * p.x + b * p.y + c) / sqrt(a * a + b * b);
		return D;
	}
	Point PointOfCrossing(const Line& A) {
		Point POC;
		POC.x = (b * A.c - A.b * c) / (a * A.b - A.a * b);
		POC.y = (A.a * c - a * A.c) / (a * A.b - A.a * b);
		return POC;
	}
	bool operator == (const Line& A) {
		return (a * A.b == b * A.a && b * A.c == c * A.b && a * A.c == c * A.a);
	}
	friend istream& operator >> (istream&, Line&);
	friend ostream& operator << (ostream&, const Line&);
	friend class Circle;
};
istream& operator >> (istream& in, Line& a) {
	in >> a.a >> a.b >> a.c;
	return in;
}
ostream& operator << (ostream& out, const Line& a) {
	out << a.a << ' ' << a.b << ' ' << a.c;
	return out;
}

class Segment {
private:
	Point a, b;
public:
	Segment() :a(0, 0), b(0, 0) {}
	Segment(const Point& p1, const Point& p2) :a(p1), b(p2) {}
	Segment(const Segment& AB) {
		a = AB.a;
		b = AB.b;
	}
	Segment& operator=(const Segment& AB) {
		a = AB.a;
		b = AB.b;
		return *this;
	}
	bool isPointOnSegment(const Point& p) {
		Line X(a, b);
		Vector AP(a, p);
		Vector BP(b, p);
		Vector AB(a, b);
		Vector BA(b, a);
		bool result = X.isPointOnLine(p) && (AP(AB) >= 0) && (BP(BA) >= 0);
		return result;
	}
	double PointDistance(const Point& p) {
		if (Vector(a, b)(Vector(a, p)) > 0 && Vector(b, a)(Vector(b, p)) > 0) {
			double D = Line(a, b).PointDistance(p);
			return D;
		}
		else if ((Vector(a, b)(Vector(a, p)) < 0 && Vector(b, a)(Vector(b, p)) > 0)) {
			double D = a.findDistance(p);
			return D;
		}
		else if ((Vector(a, b)(Vector(a, p)) > 0 && Vector(b, a)(Vector(b, p)) < 0)) {
			double D = b.findDistance(p);
			return D;
		}
		else
			return 0;
	}
	double SegmentDistance(Segment& AB) {
		double D1 = min(AB.PointDistance(a), AB.PointDistance(b));
		double D2 = min(Segment::PointDistance(AB.a), Segment::PointDistance(AB.b));
		double D = min(D1, D2);
		return D;
	}
	Point CrossingWithSegment(const Segment& AB) {
		Line X(a, b); Line Y(AB.a, AB.b);
		Point CrossPoint = X.PointOfCrossing(Y);
		if (Segment::isPointOnSegment(CrossPoint))
			return CrossPoint;
		else
			return Point(1000000, 100000);
	}
	friend istream& operator >> (istream&, Segment&);
	friend ostream& operator << (ostream&, const Segment&);
};
istream& operator >> (istream& in, Segment& AB) {
	in >> AB.a >> AB.b;
	return in;
}
ostream& operator << (ostream& out, const Segment& AB) {
	out << AB.a << ' ' << AB.b;
	return out;
}

class Circle {
private:
	double xc, yc, R;
public:
	Circle(double xc = 0, double yc = 0, double R = 0) {
		this->xc = xc;
		this->yc = yc;
		this->R = R;
	}
	Circle(const Circle& k) {
		xc = k.xc;
		yc = k.yc;
		R = k.R;
	}
	/*Circle(const vector<Point>& p) : xc(0), yc(0), R(0) { // Build circle area on known massive of point
		bool AllPointInCircleArea = true;
		for (int i = 0; i < p.size(); i++) {
			for (int j = i + 1; j < p.size(); j++) {
				for (int k = j + 1; k < p.size(); k++) {
					Circle bufer;
					Vector p1p2(p[i], p[j]), p2p3(p[j], p[k]);
					p1p2 = p1p2 / 2; p2p3 = p2p3 / 2;
					Line a1(p[i], p1p2), a2(p[j], p2p3);
					Point O = a1.PointOfCrossing(a2);
					bufer.xc = O.x; bufer.yc = O.y;
					bufer.R = sqrt((p[k].x - bufer.xc) * (p[k].x - bufer.xc) + (p[k].y - bufer.yc) * (p[k].y - bufer.yc));
					for (int l = 0; l < p.size(); l++) {
						if (bufer.PointInCircleArea(p[l]) || bufer.PointOnCircle(p[l]))
							AllPointInCircleArea *= true;
						else 
							AllPointInCircleArea *= false;
					}
					if (AllPointInCircleArea) {
						xc = bufer.xc;
						yc = bufer.yc;
						R = bufer.R;
						break;
					}
					else
						AllPointInCircleArea = true;
				}
				if (AllPointInCircleArea && xc != 0 && yc != 0 && R != 0)
					break;
			}
			if (AllPointInCircleArea && xc != 0 && yc != 0 && R != 0)
				break;
		}
	}*/
	Circle& operator=(const Circle& k) {
		xc = k.xc;
		yc = k.yc;
		R = k.R;
		return *this;
	}
	double LineDistanceFromCenter(Line& A) {
		return A.PointDistance(Point(xc, yc));
	}
	double LineDistanceFromCircle(Line& A) {
		double D = Circle::LineDistanceFromCenter(A) - R;
		return D;
	}
	pair<Point,Point> CrossingWithLine(Line& A) {
		double D = Circle::LineDistanceFromCenter(A);
		if (D > R) {
			return pair<Point,Point>(Point(1000000,1000000), Point(1000000,1000000));
		}
		else {
			Line A1(-A.b, A.a, A.b * xc - A.a * yc);
			Point O = A.PointOfCrossing(A1);
			Vector U1(O, Point(-A.b + O.x, A.a + O.y));
			Vector U2(O, Point(A.b + O.x, -A.a + O.y));
			U1 = (U1 / U1.size()) * sqrt(R * R - D * D);
			U2 = (U2 / U2.size()) * sqrt(R * R - D * D);
			return pair<Point, Point>(Point(O.x + U1.x, O.y + U1.y), Point(O.x + U2.x, O.y + U2.y));
			/*3 4 5 1 0 -8
			  2
			  11.000000 8.000000
			  0.000000 -0.000000
			  Как-то раз вывела это программы(truth)
			*/
		}
	}
	pair<Point, Point> CrossingWithCircle(Circle& k) {
		double D = Point(xc, yc).findDistance(Point(k.xc, k.yc));
		if (D == 0 && R == k.R)
			return pair<Point, Point>(Point(-1000000, -1000000), Point(-1000000, -1000000));
		else if ((D > R + k.R) || (D + k.R < R) || (D + R < k.R))
			return pair<Point, Point>(Point(1000000, 1000000), Point(1000000, 1000000));
		else if ((D == R + k.R) || (D + k.R == R) || (D + R == k.R)) {
			Vector center(Point(xc, yc), Point(k.xc, k.yc));
			center = (center / center.size()) * R;
			return pair<Point, Point>(Point(xc + center.x, yc + center.y), Point(xc + center.x, yc + center.y));
		}
		else {
			double d1 = (R * R - k.R * k.R + D * D) / (2 * D);
			double h = sqrt(R * R - d1 * d1);
			Vector center(Point(xc, yc), Point(k.xc, k.yc));
			center = (center / center.size()) * d1;
			Point O(xc + center.x, yc + center.y);
			Vector U1(-center.y, center.x);
			Vector U2(center.y, -center.x);
			U1 = (U1 / U1.size()) * h;
			U2 = (U2 / U2.size()) * h;
			return pair<Point, Point>(Point(O.x + U1.x, O.y + U1.y), Point(O.x + U2.x, O.y + U2.y));
		}
	}
	bool PointInCircleArea(const Point& p) {
		return (p.x - xc) * (p.x - xc) + (p.y - yc) * (p.y - yc) < R * R;
	}
	bool PointOnCircle(const Point& p) {
		return (p.x - xc) * (p.x - xc) + (p.y - yc) * (p.y - yc) == R * R;
	}
	pair<Line, Line> buildTangent(Point& p, pair<Point, Point>& para) {
		if (Circle::PointInCircleArea(p))
			return pair<Line, Line>(Line(1, -1, 100000), Line(1, -1, 100000));
		else if (Circle::PointOnCircle(p)) {
			double D = p.findDistance(Point(xc, yc));
			double X = sqrt(D * D - R * R);
			Circle Big(p.x, p.y, X);
			pair<Point, Point> CrossPoint = Circle::CrossingWithCircle(Big);
			para = CrossPoint;
			Line Center(Point(xc, yc), p);
			pair<Line, Line> Tangents(Line(-Center.b, Center.a, Center.b*p.x - Center.a*p.y), Line(-Center.b, Center.a, Center.b * p.x - Center.a * p.y));
			return Tangents;
		}
		else {
			double D = p.findDistance(Point(xc, yc));
			double X = sqrt(D * D - R * R);
			Circle Big(p.x, p.y, X);
			pair<Point, Point> CrossPoint = Circle::CrossingWithCircle(Big);
			para = CrossPoint;
			pair<Line, Line> Tangents(Line(p, CrossPoint.first), Line(p, CrossPoint.second));
			return Tangents;
		}
	}
	friend istream& operator >> (istream&, Circle&);
	friend ostream& operator << (ostream&, const Circle&);
};
istream& operator >> (istream& in, Circle& k) {
	in >> k.xc >> k.yc >> k.R;
	return in;
}
ostream& operator << (ostream& out, const Circle& k) {
	out << k.xc << ' ' << k.yc << ' ' << k.R;
	return out;
}

class Polygon {
private:
	vector<Point> Figure;
public:
	Polygon(int n = 0) {
		Figure.resize(n, Point(0,0));
	}
	Polygon(const vector<Point>& F2) {
		Figure.resize(F2.size());
		for (int i = 0; i < Figure.size(); i++) 
			Figure[i] = F2[i];
	}
	Polygon(const Polygon& P) {
		Figure.resize(P.Figure.size());
		for (int i = 0; i < Figure.size(); i++)
			Figure[i] = P.Figure[i];
	}
	Polygon& operator=(const Polygon& P) {
		Figure.clear();
		Figure.resize(P.Figure.size());
		for (int i = 0; i < Figure.size(); i++)
			Figure[i] = P.Figure[i];
	}
	double findSquare() {
		double S = 0;
		int i;
		for (i = 1; i < Figure.size(); i++) {
			Vector Q1(Point(0, 0), Figure[i]);
			Vector Q2(Point(0, 0), Figure[i - 1]);
			S += Q1[Q2];
		}
		Vector Q1(Point(0, 0), Figure[0]);
		Vector Q2(Point(0, 0), Figure[i - 1]);
		S += Q1[Q2];
		S *= 0.5;
		S = abs(S);
		return S;
	}
	bool isPointInPolygonArea(const Point& p) {
		double sum_corner = 0;
		int i;
		for (i = 1; i < Figure.size(); i++) {
			Vector Q1(p, Figure[i - 1]);
			Vector Q2(p, Figure[i]);
			sum_corner += atan2(Q1[Q2], Q1(Q2));
		}
		Vector Q1(p, Figure[i - 1]), Q2(p, Figure[0]);
		sum_corner += atan2(Q1[Q2], Q1(Q2));
		return sum_corner != 0;
	}
	bool isConvexPolygon() {
		bool answer = true;
		if (Vector(Figure[0], Figure[1])[Vector(Figure[1], Figure[2])] >= 0) {
			for (int i = 2; i < Figure.size(); i++) {
				answer *= (Vector(Figure[i - 2], Figure[i - 1])[Vector(Figure[i - 1], Figure[i])] >= 0);
			}
			answer *= (Vector(Figure[Figure.size() - 1], Figure[0])[Vector(Figure[0], Figure[1])] >= 0);

		}
		else {
			for (int i = 2; i < Figure.size(); i++) {
				answer *= (Vector(Figure[i - 2], Figure[i - 1])[Vector(Figure[i - 1], Figure[i])] <= 0);
			}
			answer *= (Vector(Figure[Figure.size() - 1], Figure[0])[Vector(Figure[0], Figure[1])] <= 0);

		}
		return answer;
	}
	friend istream& operator >> (istream&, Polygon&);
	friend ostream& operator << (ostream&, const Polygon&);
};
istream& operator >> (istream& in, Polygon& P) {
	for (int i = 0; i < P.Figure.size(); i++)
		in >> P.Figure[i];
	return in;
}
ostream& operator << (ostream& out, const Polygon& P) {
	cout << "Polygon P:" << endl;
	for (int i = 0; i < P.Figure.size(); i++)
		cout << "   P" << i + 1 << ": " << P.Figure[i] << endl;
	return out;
}

int main() {
	int n; cin >> n;
	Polygon P(n); cin >> P;
	bool query = P.isConvexPolygon();
	cout << (query ? "YES" : "NO");
	return 0;
}
