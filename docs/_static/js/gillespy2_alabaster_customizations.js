function setupAutodocPy() {
  const paramElements = document.querySelectorAll('.py .sig-param')
  
  Array(...paramElements).forEach((element) => {
    let brElement = document.createElement('br')
    element.parentNode.insertBefore(brElement, element)
  })
  
  const lastParamElements = document.querySelectorAll('.py em.sig-param:last-of-type')
  
  Array(...lastParamElements).forEach((element) => {
    console.log("Test")
    let brElement = document.createElement('br')    
    element.after(brElement)
  })
}

document.addEventListener('DOMContentLoaded', function() {
  console.log("Custom theme loaded.")
  
  setupAutodocPy()
})
