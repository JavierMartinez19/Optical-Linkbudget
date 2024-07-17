class MyCalendar():
    
    def __init__(self):
        pass
    
    #Leap year (año bisiesto)
    def is_leap(self, year):
        if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0): 
            leap = True
        else:
            leap = False        
        return leap
    
    def frac2date(self, frac):
        frac = str(frac)
        year = 2000 + int(frac[0:2])
        days_fractional = float(frac[2:len(frac)])
        #print("Frac: "+frac + "\nDays fractional: " + str(days_fractional))
        hours =  24 * (days_fractional - int(days_fractional))
        minutes = 60 * (hours - int(hours))
        seconds = 60 * (minutes - int(minutes))
        
        month, day = self.number2month_day(year, int(days_fractional))
        hours = int(hours)
        minutes = int(minutes)
        seconds = int(seconds)
    
        #print(str(days_fractional)+": " +str(month)+"-"+str(day)+" "+str(hours)+":"+str(minutes)+":"+str(seconds)) 
        
        return year, month, day, hours, minutes, seconds
    
    def number2month_day(self, year, day):
        if self.is_leap(year) == False:
            days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            days_in_month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]   
            
        month = 1
        while day > days_in_month[month - 1]:
            day = day - days_in_month[month - 1]
            month = month + 1  
        return month, day
    
    def float2mins(self, float):
        mins = int(float)
        secs = int((float - mins) * 60)
        return f"{mins}:{secs:02d}"

    
if __name__ == "__main__":
    mc = MyCalendar()
    year, month, day, hours, minutes, seconds = mc.frac2date(24017.65783675)
    print("{}-{}-{} {}:{}:{}".format(year, month, day, hours, minutes, seconds))