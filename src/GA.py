# Created on : Mar 27, 2024, 5:02:24 PM
# Author     : michaellevin
import sys
#print(sys.executable)

from src import Network
from src import Zone
from src import Individual
import random

class GA:
    def __init__(self, network, inflate_costs):
        self.network = network
        self.population_size = 100
        self.generations = 20
        self.mutation_rate = 0.2
        self.string_len = 6;
        
        
        self.varlinks = []
        
        for a in self.network.links:
            #print(a, a.cost)
            if a.cost > 1e-6:
                
                self.varlinks.append(a)
                
        self.g = {a:a.cost * inflate_costs for a in self.network.links}
        
        for a in self.network.links2:
            a.y = 0
            a.add_cap = 0
        
    def solve(self):
        
        print("generation", "best fitness")
        
        population = self.create_initial_population(self.population_size)
        
        for generation in range(self.generations):
            for ind in population:
                if ind.fitness > 1e10:
                    ind.fitness = self.fitness_function(ind)
                    
                    
            

            # Store the best performer of the current generation
            best_individual = self.findBest(population)
            best_fitness = best_individual.fitness
            #table.add_row([generation + 1, best_individual[0], best_individual[1], best_individual[2], best_fitness])
            print(generation+1, best_fitness)

            population = self.selection(population)

            next_population = []
            for i in range(0, len(population), 2):
                parent1 = population[i]
                parent2 = population[i + 1]

                child1, child2 = self.crossover(parent1, parent2)

                next_population.append(self.mutation(child1, self.mutation_rate))
                next_population.append(self.mutation(child2, self.mutation_rate))

            # Replace the old population with the new one, preserving the best individual
            next_population[0] = best_individual
            population = next_population
    
    def findBest(self, population):
        best = None
        best_fitness = 1e15
        
        for ind in population:
            if ind.fitness < best_fitness:
                best = ind
                best_fitness = ind.fitness
        
        
    
        
        return best
        
    def selection(self, population, tournament_size=3):
        selected = []
        for i in range(len(population)):
            winner = None
            best_fitness = 1e15
            
            for j in range(0, tournament_size):
                select = population[random.randint(0, len(population)-1)]
                
                if select.fitness < best_fitness:
                    winner = select
                    best_fitness = select.fitness
                
            selected.append(winner)
        return selected
    
    def fitness_function(self, individual):
        y = individual.calcY(self.varlinks, self.string_len)
        
        for a in self.varlinks:
            a.add_cap = y[a]
            
        self.network.tapas("UE", None)
        xhat = {a:a.x for a in self.network.links}
        obj_f = self.calcOFV(xhat, y)
        
        return obj_f
        
        
    def calcOFV(self, x, y):
        output = 0
        
        for a in self.network.links:
            y_ext = 0
            
            if a in self.varlinks:
                y_ext = y[a]
            
            output += x[a] * a.getTravelTimeC(x[a], y_ext, "UE")
            
        for a in self.varlinks:
            output += self.g[a] * y[a]
        
        return output    
        
    def create_initial_population(self, size):
        population = []
        for _ in range(size):
            individual = Individual.Individual(self.varlinks)
            individual.randomize(self.varlinks, self.string_len)
            population.append(individual)
            
        return population
                
    def crossover(self, parent1, parent2):
        
        crosspoint = random.randint(0, len(self.varlinks)-1)
        
        child1data = {a: 0 for a in self.varlinks}
        child2data = {a: 0 for a in self.varlinks}
        
        for i in range(0, crosspoint):
            a = self.varlinks[i]
            child1data[a] = parent2.data[a]
            child2data[a] = parent1.data[a]
            
        for i in range(crosspoint+1, len(self.varlinks)):
            a = self.varlinks[i]
            child1data[a] = parent1.data[a]
            child2data[a] = parent2.data[a]
            
        crosslink = self.varlinks[crosspoint]
        binary_parent1 = parent1.binary_string(crosslink, self.string_len)
        binary_parent2 = parent2.binary_string(crosslink, self.string_len)
        
        
        string_cross = random.randint(0, len(binary_parent1)-1)
        
        binary_child1 = binary_parent2[0:string_cross] + binary_parent1[string_cross:]
        binary_child2 = binary_parent1[0:string_cross] + binary_parent2[string_cross:]
        
        child1data[crosslink] = int(binary_child1, 2)
        child2data[crosslink] = int(binary_child2, 2)
        
        
        child1 = Individual.Individual(self.varlinks)
        child2 = Individual.Individual(self.varlinks)
        
        child1.setData(child1data)
        child2.setData(child2data)
        
        return child1, child2
            
    def mutation(self, individual, mutation_rate):
        for a in self.varlinks:
            if random.random() < mutation_rate:

                binary_string = individual.binary_string(a, self.string_len)
                
                mutate_point = random.randint(0, len(binary_string)-1)
                
                
                new_string = ""
                if binary_string[mutate_point] == '0':
                    new_string = binary_string[0:mutate_point]+"1"+binary_string[mutate_point+1:]  
                else:
                    new_string = binary_string[0:mutate_point]+"0"+binary_string[mutate_point+1:] 
                    
                individual.data[a] = int(new_string, 2)
                
        return individual
