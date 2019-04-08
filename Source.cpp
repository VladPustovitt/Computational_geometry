#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class Figure {
public:
	virtual double area() const = 0;
	virtual double perimeter() const = 0;
	virtual void draw() const = 0;
};

class Rectangle : public Figure {
protected:
	int h, w;
public:
	Rectangle() {
		h = 0;
		w = 0;
	}
	Rectangle(int h, int w) {
		this->h = h;
		this->w = w;
	}
	int width() const {
		return w;
	}
	int height() const {
		return h;
	}
	double area()const {
		return h * w;
	}
	double perimeter()const {
		return 2 * (h + w);
	}
	void draw() const {
		for (int y = 0; y < w; y++) {
			for (int x = 0; x < h; x++) {
				cout << '#';
			}
			cout << endl;
		}
	}
	friend istream& operator>> (istream&, Rectangle&);
	friend ostream& operator<< (ostream&, const Rectangle&);
};
istream& operator>> (istream& in, Rectangle& a) {
	in >> a.w >> a.h;
	return in;
}
ostream& operator<< (ostream& out, const Rectangle& a) {
	out << a.w << a.h;
	return out;
}

class Circle : public Figure {
protected:
	int r;
public:
	Circle() {
		r = 0;
	}
	Circle(int r) {
		this->r = r;
	}
	int radius() const {
		return r;
	}
	double area() const {
		return 3.14159265* r* r;
	}
	double perimeter() const {
		return 2 * 3.14159265* r;
	}
	void draw() const {
		for (int x = 0; x < 2 * r + 1; x++) {
			for (int y = 0; y < 2 * r + 1; y++) {
				cout << ((x - r) * (x - r) + (y - r) * (y - r) <= r * r ? '#' : '.');
			}
			cout << endl;
		}
	}
	friend istream& operator>> (istream&, Circle&);
	friend ostream& operator<< (ostream&, const Circle&);
};
istream& operator>> (istream& in, Circle& a) {
	in >> a.r;
	return in;
}
ostream& operator<< (ostream& out, const Circle& a) {
	out << a.r;
	return out;
}

class HollowRectangle : public Rectangle {
public:
	HollowRectangle(int h, int w) : Rectangle(h, w) {}
	void draw() const {
		for (int x = 0; x < w; x++) {
			for (int y = 0; y < h; y++) {
				cout << (x == 0 || x == w - 1 || y == 0 || y == h - 1 ? '#' : '.');
			}
			cout << endl;
		}
		
	}
};

class HollowCircle : public Circle {
public:
	HollowCircle(int r) : Circle(r) {}
	bool p(int a, int b) const {
		return (((a - r) * (a - r) + (b - r) * (b - r)) <= r * r );
	}
	void draw() const {
		for (int y = 0; y < 2 * r + 1; y++) {
			for (int x =  0; x < 2 * r + 1; x++) {
				cout << (p(x,y) && (!p(x-1,y) || !p(x,y-1) || !p(x+1,y) || !p(x, y+1)) ? '#' : '.');
			}
			cout << endl;
		}
	}
};

double array_area(Figure** a, int n) {
	double sum_area = 0;
	for (int i = 0; i < n; i++) {
		sum_area += a[i]->area();
	}
	return sum_area;
}
double array_perimeter(Figure** a, int n) {
	double sum_perimeter = 0;
	for (int i = 0; i < n; i++) {
		sum_perimeter += a[i]->perimeter();
	}
	return sum_perimeter;
}

void draw_array(Figure** a, int n) {
	for (int i = 0; i < n; i++) {
		a[i]->draw();
		cout << endl;
	}
}

int count_Rectangle(Figure** a, int n) {
	int sum = 0;
	for (int i = 0; i < n; i++) {
		Rectangle* p;
		p = dynamic_cast<Rectangle*>(a[i]);
		if (p != NULL)
			sum++;
	}
	return sum;
}

int count_Circle(Figure** a, int n) {
	int sum = 0;
	for (int i = 0; i < n; i++) {
		Circle* p;
		p = dynamic_cast<Circle*>(a[i]);
		if (p != NULL)
			sum++;
	}
	return sum;
}

int main() {
	int n;
	cin >> n;
	Figure** fig = new Figure* [n];
	string t;
	for (int i = 0; i < n; i++) {
		cin >> t;
		if (t == "Rectangle") {
			int w, h;
			cin >> w >> h;
			fig[i] = new Rectangle(w, h);
		}
		else if (t == "Circle") {
			int r;
			cin >> r;
			fig[i] = new Circle(r);
		}
		else if (t == "HollowRectangle") {
			int w, h;
			cin >> w >> h;
			fig[i] = new HollowRectangle(w, h);
		}
		else if (t == "HollowCircle") {
			int r;
			cin >> r;
			fig[i] = new HollowCircle(r);
		}
	}
	draw_array(fig, n);
	return 0;
}