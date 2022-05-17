const canvas = document.getElementById('canvas')
canvas.style.width = canvas.width / 2 + 'px'
canvas.style.height = canvas.height / 2 + 'px'

const kernel = new GPU2D(canvas)

for (const input of document.getElementsByTagName('input')) {
  let uniformName = input.id.split('-').join('.')
  uniformName = uniformName.replace('planetcenter', 'planetCenter')
  input.addEventListener('input', e => {
    console.log(uniformName, input.value)

    if (uniformName.endsWith('.x') || uniformName.endsWith('.y') || uniformName.endsWith('.z')) {
      let vecUniformName = uniformName.slice(0, -2)

      let arr = [0, 0, 0]
      let vecId = input.id.slice(0, -2)
      for (let i = 0; i < 3; i++) {
        let id = vecId + '-' + 'xyz'[i]
        let input = document.getElementById(id)
        arr[i] = parseFloat(input.value)
      }

      kernel.setUniformValue(vecUniformName, (gl, l) => gl.uniform3fv(l, arr))
    } else {
      kernel.setUniformValue(uniformName, (gl, l) => gl.uniform1f(l, +e.target.value))
    }
  })
}

fetch('atmosphere.vert').then(response => response.text()).then(code => {
  const program = kernel.create2dProgram(code)
  kernel.useProgram(program)

  kernel.setUniformValue('sun.center', (gl, l) => gl.uniform3f(l, 0, 1.5, -3))
  kernel.setUniformValue('sun.radius', (gl, l) => gl.uniform1f(l, 0.5))
  kernel.setUniformValue('planetCenter', (gl, l) => gl.uniform3f(l, 0, -0.801, 0))

  requestAnimationFrame(render)
})

let prev = 0
function render(now) {
  let delta = now - prev
  prev = now
  const fps = 1000 / delta
  document.getElementById('fps').innerText = fps.toFixed(2)

  kernel.draw()
  
  requestAnimationFrame(render)
}