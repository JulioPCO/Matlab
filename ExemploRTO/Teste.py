class levedura :
    def __init__ (self,name,mass):
        self.name = name
        self.mass = mass

class saccharomyces(levedura):
    def __init__(self,name,mass,glycose,rate):
        super().__init__(name,mass)
        self.glycose = glycose
        self.rate = rate

class cerevisiae(saccharomyces):
    food = True
    def __init__(self,name,mass,glycose,rate,efficacy):
        super().__init__(name,mass,glycose,rate)
        self.efficacy = efficacy

    def batch_reproduction_rate(self):
        self.reprod_rate = self.mass * self.glycose * self.efficacy/100000
        return self.reprod_rate 

    def cstr_reprod_rate(self):
        pass


class bioma:
    def __init__(self,name,max_pop):
        self.name = name
        self.max_pop = max_pop
        self.bioma_levedura = []

    def add_levedura(self, leved):
        if len(self.bioma_levedura) < self.max_pop:
            self.bioma_levedura.append(leved)
        else:
            print('Full')

levedura1 = levedura('lv1',30)

print(levedura1.name)

saccharomyces1 = saccharomyces('sc1',40,1000,0.47)

print(saccharomyces1.rate)

saccerevisiae1 = cerevisiae('cr1',40,1000,0.47,0.8)

print(saccerevisiae1.efficacy)

if saccerevisiae1.food == True : 
    print(saccerevisiae1.batch_reproduction_rate())

bioma1 = bioma('bio1',3)

bioma1.add_levedura(saccerevisiae1)

print(bioma1.bioma_levedura)

print(bioma1.bioma_levedura[0].name)

objects = []
for i in range(2,5):
    objects.append(cerevisiae('cr1',40+10*i,1000,0.47,0.8))

print(objects)

for i in objects:
    bioma1.add_levedura(i)

print(bioma1.bioma_levedura)    