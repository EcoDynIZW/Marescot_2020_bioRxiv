;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;------------------------------------INITIALISATION------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

breed [hyenas hyena]

patches-own
[
  clan-territory ;; true if the patch hold a clan territory false otherwise
  infeclan   ;; true if there is at least one infected agent in the clan false otherwise
  counterinfec
]

globals
[
 abundance ; mean abundance over the disease course and abundance at the end of a time step (week) when all biological processes occured
 infecpatch ; indicated whether a patch is infected or not
 clan-number ; determined from the carrying capacity abd the maximum group size.
 clan-ext ;
 weeklysurv ; convert the yearly survival to weekly survival

 ;----------initialization---------------------

 meanIL ; effective infection length
 infectivity ; infectivity of the pathogen adjusted to control for R0
 young ;proportion of juveniles staying in comunal den

 ;----------Output globals---------------------
 adults ;number of adults
 timetoextinction ; record the week at which the virus went extinct.
 ext; is a transitory variable so timetoextinction does not keep increasing over time it also indicates whether or not the virus went extinct
 minabundance ; minimum abundane over the disease course
 maxprev inf ;maximum prevalence during the disease course. "inf" is just a temporary variable that is updated every time prevalence overpasses the maximum value
 sealed ; number of sealed group in which all adults are recovered
 open; number of unlocked groups
 locked; number of in
 aveimeff immefficiency; average proportion of susceptible individuals in the dens in sealed groups. "immefficiency is just a temporary variable that is updated every time efficiency overpasses the maximum value
 nonprot ; proportion of susceptible in the population that are not protected
 infecub
 cub
 cuboutput
 tickyear
]

hyenas-own
[
 age                 ;; age of individuals in weeks
 clan.ID             ;; Identity of the clan (pacth gathering individuals of a given clan)
 disease             ;; "Susceptible" "infected", "recovered"
 virus-check-timer   ;; number of ticks since an individual first caught the virus
 become-infected     ;; variable used to updated the infection of individuals which has to occure synchronously at each weekly step
                     ;; otherwise the infection probability which is a function of the number of infected individuals varies between individuals (depend on the order of individual loop)
 initial-infected    ;; record the indivduals that are the first to be infected when a pathogen invade
]

 ;;extensions[profiler]




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;------------------------------------SETUPT------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup ; initialization-------------------------------------------------------------------

  clear-all
  ; random-seed 80
  set clan-number carrying-capacity / group-size ; calculate the number of patches

  ;--------- Distribution of clans across the grid landscape----------
  ask n-of clan-number patches
  [
    ; each clan has an average size of average_clan
    sprout-hyenas (group-size)
    [
      set age random-exponential (104)  ; randomly assign an age to individuals from 0 to 10 years old
      set disease "susceptible" ; set the entire clans to be susceptible
      set virus-check-timer 0
      set color green ; green to represent susceptible individuals
      set pcolor green ; green to represent susceptible patches
      ask patch-here
      [
        set clan-territory true ; store the information that this pacth has a settled clan
        set infeclan false ; making sure that all group are healthy at start
      ]
    ]
  ]

  ; Ask one individual to be infected by the virus. this indivudal must me part of the immunity barrier and not young individuals in the dens

  ask one-of hyenas
  [
    set disease "infected"
    set virus-check-timer 1 ; first week of infection
    set color red           ; red for infected individuals
  ]


  ; initialize globals
  set minabundance carrying-capacity ; initialize minimum abundance
  set maxprev 0 ; start simulation with a single infected individuals
  set infecpatch 1 ; one patch is infected
  set ext false ; make sure the check variable for virus extinction is false
  set timetoextinction 1000



  ;---------------seting for calculation of mean number ofoffspring per breeding adult

   set cub 0
   set cuboutput 0
   set tickyear 0
   
;---------------calculation of infectivity given the Rnot

  set weeklysurv   survival ^ (1 / 52)

  let listyoung n-values (immunity-barrier + 1) [ ?1 -> ?1 + 1 ]
  foreach listyoung
  [ ?1 ->
    set young young + (1  - weeklysurv) * weeklysurv ^ (?1 - 1)
  ]

  let listval n-values (infection-length) [ ?1 -> ?1 ]

  foreach listval
  [ ?1 ->
    set meanIL  meanIL + ((1 - ((1 - virulence) * weeklysurv)) * (weeklysurv * (1 - virulence))^ ?1 ) * ?1
  ]

  set meanIL  meanIL +  ((weeklysurv *(1 - virulence))^ infection-length) * infection-length

ifelse frequency_transmission = true ; infection process under frequency dependent transmission
  [
    set infectivity rnot / (meanIL * (within-group.contact + between-group.contact * (1 - young)))
  ]
  [ ; otherwise infection process under density dependent transmission
    set infectivity rnot / (meanIL * (within-group.contact * (group-size - 1) + between-group.contact * (carrying-capacity - group-size) * (1 - young) ^ 2))
  ]

  reset-ticks ; reinitialize the model when clicking on set-up

end ;----------------------------------------------------------------------------------------




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;------------------------------------PROCESSES------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; run the successive biological process we defined in seperate function
; We assumed processes to occur as follows, infection, birth, survival, dispersal and growth

;----------------------------------RUN-------------------------------------------------------------
to go
  ;;profiler:start

  infection ; infection process
  do-virus-checks ; procedure which track the time since infection and assign an immune status to individuals who survived infection until the maximum infection-length
  death ; natural death and disease-related mortaltiy
  reproduce ; reproduction

  ifelse fission = true
  [
     budding
  ]
  [
     migration
  ]
  aging

  outputvar ; store the response variable we are interested in. abundance and maximum prevalence
  ;getimmunity ; calculate whether an immunity barrier emerge and how efficient it can be

  tick
   ;;profiler:stop
   ;;print profiler:report
   ;;profiler:reset

end




;;------------------------------------INFECTION------------------------------------------


to infection ; modeling the infecion process

  let popsize count hyenas ; total abundance at the moment when infection occurs

  ask patches with [clan-territory = true] ; select patches occupied by a clan
  [
    let clansize count other hyenas-here ; total number of individuals in the group at the moment of infection
    ifelse  clansize > 0
    [
      ;;;;;;; to calculate the between group per capita contact rate (frequency dependent so a reduced between group contact does not depent on age-structured interactions)

      ; total of individual outside the group by counting all individuals located at a different distance than "myseld" as all individuals within the groups have the same coordinates
      let otherclan-pop-size count hyenas with [distance myself != 0 and age >= immunity-barrier]
      if  otherclan-pop-size = 0 [set  otherclan-pop-size 1] ; to avoid message error due to division by 0 when calculating the frequency dependent infection
      ; count the number of infected adults outside the clan that can potenitally infect the hyneas of clan
      let without-hyenas count hyenas with [disease = "infected" and age >= immunity-barrier and distance myself != 0]

      ;;;;;;;; to calculate the within group per capita contact rate (frequency dependent so within group contact only does not depend on group size)

      let within-hyenas count hyenas-here  with [disease = "infected"]  ; count the number of infected individuals in the group that can potentially infect susceptibles
      ifelse within-hyenas > 0 [set infeclan true] [set infeclan false] ; checking the disease status of the clan
                                                                        ; uncomment if we want the between transmisison to be a function of a distance

      ;;;;;;;;;;;;;;;;;;;;;;;;Infection of Independent;;;;;;;;;;;;;;;;;;;;

      ;within and between group infection of adults in the den
         ask hyenas-here with [disease = "susceptible" and age >= immunity-barrier] ; loop over suceptible individuals of the clan that are at the immunity barrier and can be infected from outside
      [
       ifelse frequency_transmission = true ; infection process under frequency dependent transmission
       [
        if random-float 1 < (1 - 1 / exp((infectivity * within-group.contact * within-hyenas / clansize) +  (infectivity * between-group.contact * without-hyenas  / otherclan-pop-size)))
        [
          set become-infected true ; if the individual was infected by either an individual from the same clan or from outside
        ]
       ]
       [
        if random-float 1 < (1 - 1 / exp((infectivity * within-group.contact * within-hyenas) +  (infectivity * between-group.contact * without-hyenas)))
        [
          set become-infected true ; if the individual was infected by either an individual from the same clan or from outside
        ]
       ]

      ]

      ;;;;;;;;;;;;;;;;;;;;;;;;Infection of Young;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;only within-group infection of youngs in the dens

      ask hyenas-here with [disease = "susceptible" and age < immunity-barrier] ; model the infection of individual within the communal den that can only be infected by members of the clan
      [
        ifelse frequency_transmission = true ; infection process under frequency dependent transmission
       [

        if random-float 1 < (1 - 1 / exp((infectivity * within-group.contact * within-hyenas / clansize)))
        [
          set become-infected true
        ]
       ]
       [
        if random-float 1 < (1 - 1 / exp((infectivity * within-group.contact * within-hyenas)))
        [
          set become-infected true
        ]
       ]
      ]
    ]

    [ ; if clan size = 0  setting the patch variable to extinct if the group went extinct

      set clan-ext  "extinct"
    ]

  ]

  ; Synchronization of the infection of individuals at each weekly step
  ; otherwise the infection probability depend on the order of individual being processed

    ask hyenas with [become-infected = true]
    [
      set disease "infected"
      set become-infected  false
      set color red
    ]

end


to do-virus-checks ; function that change individuals disease status synchronously in CDV individuals gain life immunity after 4 weeks of infection periods

  ask hyenas with [disease = "infected" and virus-check-timer = infection-length]
  [
    set disease "recovered"
    set color grey
  ]

  ask hyenas with [disease = "infected"]
  [
    set virus-check-timer virus-check-timer + 1
  ]

end


;;------------------------------------SURVIVAL------------------------------------------

to death
  ask hyenas
  [
    if(random-float 1) < (1 - weeklysurv) [die]

  ]

  ask hyenas with [disease = "infected"]
  [
     if random-float 1 < virulence [die]
  ]

  if (count hyenas with[disease = "infected"] = 0 and ext = false) ; check the time when pathogen goes extinct
  [
     set ext true
  ]


end

;;------------------------------------ GROWTH------------------------------------------

to aging
  ask hyenas
  [
     set age age + 1
  ]
end
;;------------------------------------REPRODUCTION------------------------------------------

to reproduce
    set tickyear tickyear + 1
    ask patches with [clan-territory = true and count hyenas-here <=  group-size]
    [
       ask hyenas-here with [age >= age-repro]
       [
        if (random-float 1) < (fertility)
        [

            set cub cub + 1
        ]
      ]
     ]

      if tickyear = 52
      [
       set cuboutput cub
       ask patches with [clan-territory = true and count hyenas-here <=  group-size]
       [
        ask n-of (count hyenas-here with [age >= age-repro] / 2) hyenas-here with [age >= age-repro]
        [
        if cuboutput > 0
        [
          hatch-hyenas 1 ;round (52 * fertility)
          [
            ;set birthdate ticks
            set age 0
            set disease "susceptible"
            set virus-check-timer 0
          ]
          set cuboutput cuboutput - 1 
        ]
        ]
       ]
        set cub 0
        set tickyear 0
      ]



end




;;------------------------------------DISPERSAL---------------------------------------------

to migration
  ;;ask patches with [clan-territory = true and count hyenas-here >  average-clan]
  ask patches with [clan-territory = true]
  [
      ask hyenas-here with [age = age-repro]
      [
          move-to one-of other patches with [clan-territory = true]
      ]

  ]
end

to budding
  ;;ask patches with [clan-territory = true and count hyenas-here >  average-clan]
  ask patches with [clan-ext = "extinct"]
  [
    let pxID pxcor
    let pyID pycor
    if(count patches with [clan-territory = true and count hyenas-here >=  group-size] > 0)
    [
      ask one-of patches with [clan-territory = true and count hyenas-here >=  group-size]
      [
        let halfsize round count hyenas-here / 2
        ask n-of halfsize  hyenas-here
        [
          move-to patch pxID pyID
        ]
      ]
    ]
    set clan-ext "occupied"
  ]
end



;;------------------------------------SUPLEMENTARY FUNCTIONS------------------------------------------

to getimmunity
 set immefficiency 0
 set infecub 0
 set sealed 0
 set open 0
 set nonprot 0


  ; variable in calcualted in sealed groups
 ask patches with [clan-territory = true and count hyenas-here with [age >= immunity-barrier] =  count hyenas-here with [age >= immunity-barrier and disease = "recovered"]]
 [
    set sealed sealed + 1 ; counting the number of sealed group
   if count hyenas-here with [age < immunity-barrier] > 0
   [

     set immefficiency immefficiency + count hyenas-here with [age < immunity-barrier and disease = "susceptible"] / count hyenas-here with [age < immunity-barrier]
     set locked locked + count hyenas-here with [age < immunity-barrier and disease != "susceptible"] / count hyenas-here with [age < immunity-barrier]

     set infecub infecub + count hyenas-here with [age < immunity-barrier and disease = "infected"] ; counting the number of infected cubs in sealed groups
   ]
 ]


 if count hyenas with [age < immunity-barrier and disease = "infected"] > 0
 [
    set infecub infecub / count hyenas with [age < immunity-barrier and disease = "infected"] ; proportion of infected cubs that are in sealed groups
 ]

 ; variable obtained in non protected susceptile in groups which are not sealed
 ask patches with [clan-territory = true and count hyenas-here with [age >= immunity-barrier] =  count hyenas-here with [age >= immunity-barrier and disease != "recovered"]]
 [
    set open open + 1
   if count hyenas-here with [age < immunity-barrier] > 0
   [

     set nonprot nonprot + count hyenas-here with [age < immunity-barrier and disease = "susceptible"]  ; proprotion of susceptibles juveniles that are not protected
   ]
 ]


 if count hyenas with [age < immunity-barrier and disease = "susceptible"] > 0
 [
   set nonprot nonprot / count hyenas with [age < immunity-barrier and disease = "susceptible"] ; proportion of susceptible juveniles that are not protected

 ]


 if sealed > 0
 [
   set immefficiency (immefficiency) / sealed ; efficiency of the immunity barrier mean proportion of susceptible protected by sealed groups
   set locked (locked) / sealed ; mean proportion of juveniles in sealed groups that are either infected or immune

  ]
 if immefficiency > aveimeff
 [
 ;  set aveimeff immefficiency ; in case we want the myamimum efficienccy over the 110 years
 ]
end

 to colorset
 ask patches with [clan-territory = true]
 [
   ifelse any? hyenas-here with [disease ="infected"]
   [
     set pcolor red
   ]
   [
     ifelse any? hyenas-here with [disease ="susceptible"]
     [
       set pcolor green
     ]
     [
        set pcolor grey
     ]
   ]
 ]
end

to outputvar
   set abundance  count hyenas
   set adults count hyenas with [age >= immunity-barrier]
   if abundance < minabundance
   [
     set minabundance abundance ; record the mimnimum abundance over the course of the disease to measure pathogen impact
   ]

   if abundance = 0
   [
     stop ; stop the simulation if the host goes extinct
   ]

   if abundance > 0
   [
     set inf count hyenas with [disease = "infected"] / abundance
     if inf > maxprev
     [
       set maxprev inf ; measures the maximum prevalence over the course of the disease
     ]
   ]

  if ext = true ; if it was extinction for the first time
  [
       set timetoextinction ticks ; uncommented if we want to record the time at which the pathogen went extinct
       set ext "extinct" ; Extintinction status of the pathoggen changing the ext variable so the timetoextinction varibale does not keep counting overtime
       stop ; stop the simulation when pathogen goes extinct
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
246
14
632
401
-1
-1
11.455
1
10
1
1
1
0
1
1
1
-16
16
-16
16
1
1
1
ticks
30.0

BUTTON
58
67
121
100
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
57
116
120
149
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
31
275
203
308
group-size
group-size
0
200
10.0
1
1
NIL
HORIZONTAL

SLIDER
32
319
204
352
carrying-capacity
carrying-capacity
0
5000
1000.0
100
1
NIL
HORIZONTAL

SLIDER
30
367
202
400
fertility
fertility
0
1
0.01
0.01
1
NIL
HORIZONTAL

SLIDER
30
418
205
451
survival
survival
0
1
0.9
0.01
1
NIL
HORIZONTAL

SLIDER
821
32
1005
65
within-group.contact
within-group.contact
0
1
1.0
0.0001
1
NIL
HORIZONTAL

SLIDER
822
72
1001
105
between-group.contact
between-group.contact
0
1
0.01
0.1
1
NIL
HORIZONTAL

SLIDER
822
117
994
150
virulence
virulence
0
1
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
823
160
995
193
immunity-barrier
immunity-barrier
0
1000
52.0
1
1
NIL
HORIZONTAL

SLIDER
825
263
997
296
infection-length
infection-length
0
100
30.0
1
1
NIL
HORIZONTAL

SLIDER
824
215
996
248
age-repro
age-repro
0
208
104.0
1
1
NIL
HORIZONTAL

SWITCH
37
220
140
253
fission
fission
1
1
-1000

SLIDER
39
165
211
198
rnot
rnot
0
50
2.0
0.1
1
NIL
HORIZONTAL

SWITCH
765
416
952
449
frequency_transmission
frequency_transmission
0
1
-1000

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment_countingoffspring_BIRTHPULSE" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="260"/>
    <metric>maxprev</metric>
    <metric>ext</metric>
    <metric>timetoextinction</metric>
    <metric>abundance</metric>
    <metric>minabundance</metric>
    <metric>infectivity</metric>
    <metric>cub</metric>
    <metric>breedingadult</metric>
    <metric>count hyenas with [disease = "infected"]</metric>
    <enumeratedValueSet variable="within-group.contact">
      <value value="0.1"/>
      <value value="0.5"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infection-length">
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="group-size">
      <value value="10"/>
      <value value="50"/>
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rnot">
      <value value="0.5"/>
      <value value="1"/>
      <value value="2"/>
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="virulence">
      <value value="0"/>
      <value value="0.05"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-barrier">
      <value value="0"/>
      <value value="52"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="clan-number">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="between-group.contact">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="survival">
      <value value="0.6"/>
      <value value="0.75"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fertility">
      <value value="0.01"/>
      <value value="0.025"/>
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fission">
      <value value="false"/>
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
