favorite_language = {
    'jen':['python','ruby'],
    'sarah':['c'],
    'edward':['ruby','go'],
    'phil':['python','haskell'],
}

for name, languages in favorite_language.items() :
    print("\n" + name.title() + "'s favorite language are:")
    for language in languages :
        print("\t" + language.title())