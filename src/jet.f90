   function jet(x, yn, h)
   implicit real*8 (a-z)
      t1 = cos(x)
      t2 = cos(yn)
      t3 = t2+t1
      t4 = sin(x)
      t5 = sin(yn)
      t6 = (-t3*t5)-t4
      t7 = t3**2
      t8 = (-t5*t6)-t2*t7-t1
      t9 = (-t5*t8)-3*t2*t3*t6+t3**3*t5+t4
      Y1 = t3
      Y2 = t6
      Y3 = t8
      Y4 = t9
      Y5 = (-3*t2*t6**2)-t5*t9-4*t2*t3*t8+6*t7*t5*t6+t2*t3**4+t1
      jet = yn+h*((h*((h*((h*(Y5*h/5.0d+0+Y4))/4.0d+0+Y3))/3.0d+0+Y2))/2.0d+0+Y1)
    end function jet
