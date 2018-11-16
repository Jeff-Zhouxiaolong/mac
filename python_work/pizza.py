pizza = {
    'crust':'thick',
    'toppings':['mushrooms','extra cheese'],
}


print("You order a " + pizza['crust'] + "-curst pizza" + "with the following toppings")

for topping in pizza['toppings'] :
    print("\t" + topping)